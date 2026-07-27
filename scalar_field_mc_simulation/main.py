import numpy as np
import matplotlib.pyplot as plt
from scalar_field import ScalarField2D

def run_simulation():
    # パラメータ設定
    L = 20          # 格子サイズ
    m = 0.5         # 質量
    delta = 1.5     # 更新ステップサイズ
    n_thermal = 500 # 熱平衡化のためのスイープ数
    n_measure = 1000 # 測定のためのスイープ数
    interval = 10   # 測定間隔

    print(f"Starting Simulation: L={L}, m={m}, delta={delta}")
    model = ScalarField2D(L, m, delta)

    # 1. 熱平衡化 (Thermalization)
    print("Thermalizing...")
    actions_thermal = []
    for s in range(n_thermal):
        acc_rate = model.metropolis_sweep()
        if s % 50 == 0:
            action = model.calculate_total_action()
            actions_thermal.append(action)
            print(f"  Step {s:4d}: Action = {action:10.2f}, Acc Rate = {acc_rate:.2f}")

    # 2. 測定 (Measurement)
    print("\nMeasuring...")
    phi2_list = []
    max_r = L // 2
    g_r_accumulator = np.zeros(max_r)
    n_samples = 0

    for s in range(n_measure):
        acc_rate = model.metropolis_sweep()
        if s % interval == 0:
            phi2 = model.measure_phi2()
            phi2_list.append(phi2)
            
            # 相関関数の計算
            g_r_accumulator += model.calculate_correlation_function(max_r)
            n_samples += 1

            if s % 100 == 0:
                print(f"  Step {s:4d}: <phi^2> = {phi2:.4f}")

    # 統計計算
    phi2_mean = np.mean(phi2_list)
    g_r_final = g_r_accumulator / n_samples

    # 相関関数から質量 m を推定する
    # G(r) ~ exp(-m_eff * r) なので、log(G(r)) ~ -m_eff * r
    # 最初の数点は格子効果が強く、後半はノイズが大きいため、中間の範囲でフィットを行う
    fit_range = slice(1, L // 4 + 1) # 例: r=1 から L/4 まで
    r_vals = np.arange(max_r)[fit_range]
    log_g = np.log(g_r_final[fit_range])
    
    # 最小二乗法で傾きを求める
    slope, intercept = np.polyfit(r_vals, log_g, 1)
    m_eff = -slope

    print(f"\nResults:")
    print(f"  Final <phi^2> = {phi2_mean:.6f}")
    print(f"  Input mass m  = {m:.4f}")
    print(f"  Estimated effective mass m_eff = {m_eff:.4f}")

    # 可視化
    plt.figure(figsize=(15, 5))

    # 1. フィールド構成
    plt.subplot(1, 3, 1)
    plt.imshow(model.phi, cmap='RdBu')
    plt.colorbar()
    plt.title("Field Configuration $\phi(x)$")

    # 2. 相関関数 G(r)
    plt.subplot(1, 3, 2)
    plt.plot(range(max_r), g_r_final, 'o-')
    plt.axhline(0, color='black', linestyle='--', alpha=0.3)
    plt.title("Correlation Function $G(r)$")
    plt.xlabel("Distance $r$")
    plt.ylabel(r"$\langle \phi(0)\phi(r) \rangle$")
    plt.grid(True, alpha=0.3)

    # 3. 相関関数の対数プロット (指数減衰の確認)
    plt.subplot(1, 3, 3)
    plt.semilogy(range(max_r), g_r_final, 'o-')
    plt.title("$G(r)$ (Log scale)")
    plt.xlabel("Distance $r$")
    plt.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    plt.savefig("correlation_results.png")
    print("\nSaved new visualization to 'correlation_results.png'")

if __name__ == "__main__":
    run_simulation()
