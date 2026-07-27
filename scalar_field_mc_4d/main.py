import numpy as np
import matplotlib.pyplot as plt
from scalar_field_4d import ScalarField4D

def run_simulation():
    # 4次元なので格子サイズを抑える (6^4 = 1296点)
    L = 6
    m = 1.0
    delta = 1.5
    n_thermal = 200
    n_measure = 400
    interval = 5

    print(f"Starting 4D Simulation (Euclidean): L={L}, m={m}")
    model = ScalarField4D(L, m, delta)

    print("Thermalizing...")
    for s in range(n_thermal):
        acc = model.metropolis_sweep()
        if s % 50 == 0:
            print(f"  Step {s:4d}: Acc Rate = {acc:.2f}")

    print("\nMeasuring...")
    max_r = L // 2
    g_r_acc = np.zeros(max_r)
    n_samples = 0

    for s in range(n_measure):
        model.metropolis_sweep()
        if s % interval == 0:
            g_r_acc += model.calculate_correlation(max_r)
            n_samples += 1

    g_r = g_r_acc / n_samples
    r_vals = np.arange(1, max_r)
    
    # 4次元自由場: G(r) ~ 1/r^2 * exp(-mr)
    # フィッティング用: log(G(r) * r^2) ~ -m*r
    gr_r2 = g_r[1:] * (r_vals**2)
    slope, _ = np.polyfit(r_vals, np.log(gr_r2), 1)
    m_eff = -slope

    print(f"\nResults:")
    print(f"  Input m = {m:.2f}, Estimated m_eff = {m_eff:.4f}")

    # 可視化
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(range(max_r), g_r, 'o-')
    plt.title("Correlation Function $G(r)$ (4D)")
    plt.xlabel("Distance $r$")
    plt.grid(True, alpha=0.3)

    plt.subplot(1, 2, 2)
    plt.semilogy(r_vals, gr_r2, 's-', label='$G(r) \cdot r^2$')
    plt.title("Log Plot: $\log(G(r) \cdot r^2)$")
    plt.xlabel("Distance $r$")
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    plt.savefig("results_4d.png")
    print(f"Saved results to Desktop/scalar_field_mc_4d/results_4d.png")

if __name__ == "__main__":
    run_simulation()
