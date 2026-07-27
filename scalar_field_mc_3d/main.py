import numpy as np
import matplotlib.pyplot as plt
from scalar_field_3d import ScalarField3D

def run_simulation():
    # パラメータ設定 (3次元なのでLを少し小さく設定)
    L = 10          # 格子サイズ (10x10x10)
    m = 0.5         # 質量
    delta = 1.2     # 更新ステップサイズ
    n_thermal = 300 # 熱平衡化のためのスイープ数
    n_measure = 500 # 測定のためのスイープ数
    interval = 5    # 測定間隔

    print(f"Starting 3D Simulation: L={L}, m={m}, delta={delta}")
    model = ScalarField3D(L, m, delta)

    # 1. 熱平衡化 (Thermalization)
    print("Thermalizing...")
    for s in range(n_thermal):
        acc_rate = model.metropolis_sweep()
        if s % 50 == 0:
            action = model.calculate_total_action()
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
    # 3次元の場合、G(r) ~ (1/r) * exp(-m_eff * r) の形になるが、
    # 簡易的に指数減衰として推定を試みる
    fit_range = slice(1, max_r)
    r_vals = np.arange(max_r)[fit_range]
    # G(r) * r ~ exp(-m_eff * r) と仮定
    log_gr = np.log(g_r_final[fit_range] * r_vals)
    
    # 最小二乗法
    slope, intercept = np.polyfit(r_vals, log_gr, 1)
    m_eff = -slope

    print(f"\nResults:")
    print(f"  Final <phi^2> = {phi2_mean:.6f}")
    print(f"  Input mass m  = {m:.4f}")
    print(f"  Estimated effective mass m_eff = {m_eff:.4f}")

    # 可視化
    plt.figure(figsize=(15, 5))

    # 1. フィールド構成 (中央の断面を表示)
    plt.subplot(1, 3, 1)
    plt.imshow(model.phi[L//2, :, :], cmap='RdBu')
    plt.colorbar()
    plt.title(f"Field Section $\phi(L/2, y, z)$")

    # 2. 相関関数 G(r)
    plt.subplot(1, 3, 2)
    plt.plot(range(max_r), g_r_final, 'o-')
    plt.axhline(0, color='black', linestyle='--', alpha=0.3)
    plt.title("Correlation Function $G(r)$")
    plt.xlabel("Distance $r$")
    plt.ylabel(r"$\langle \phi(0)\phi(r) \rangle$")
    plt.grid(True, alpha=0.3)

    # 3. 相関関数の対数プロット (G(r)*r)
    plt.subplot(1, 3, 3)
    plt.semilogy(r_vals, g_r_final[fit_range] * r_vals, 'o-', label='$G(r) \cdot r$')
    plt.title("Log plot of $G(r) \cdot r$")
    plt.xlabel("Distance $r$")
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    plt.savefig("correlation_results_3d.png")
    print("\nSaved new visualization to 'correlation_results_3d.png'")

if __name__ == "__main__":
    run_simulation()
