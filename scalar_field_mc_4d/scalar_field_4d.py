import numpy as np

class ScalarField4D:
    def __init__(self, L, m, delta=0.5):
        """
        4次元実スカラー場の格子シミュレーション (Euclidean)
        Minkowski空間の物理は、ウィック回転後のEuclidean空間で計算するのが標準的です。
        
        Args:
            L (int): 各次元の格子サイズ (L^4)
            m (float): 質量
            delta (float): メトロポリス更新の幅
        """
        self.L = L
        self.m = m
        self.delta = delta
        # 4次元配列 (x, y, z, t)
        self.phi = np.zeros((L, L, L, L))
        
    def get_local_action(self, x, y, z, t):
        """指定された点 (x, y, z, t) の局所作用を計算"""
        phi_curr = self.phi[x, y, z, t]
        
        # 8つの隣接点 (4次元なので各軸±方向)
        phi_neighbors = (
            self.phi[(x + 1) % self.L, y, z, t] + self.phi[(x - 1) % self.L, y, z, t] +
            self.phi[x, (y + 1) % self.L, z, t] + self.phi[x, (y - 1) % self.L, z, t] +
            self.phi[x, y, (z + 1) % self.L, t] + self.phi[x, y, (z - 1) % self.L, t] +
            self.phi[x, y, z, (t + 1) % self.L] + self.phi[x, y, z, (t - 1) % self.L]
        )
        
        # 4次元なので運動項の係数は 4.0
        kinetic = 4.0 * phi_curr**2 - phi_curr * phi_neighbors
        mass_term = 0.5 * (self.m**2) * (phi_curr**2)
        
        return kinetic + mass_term

    def metropolis_sweep(self):
        """1スイープ実行"""
        accepted = 0
        # 4重ループ
        for x in range(self.L):
            for y in range(self.L):
                for z in range(self.L):
                    for t in range(self.L):
                        old_val = self.phi[x, y, z, t]
                        old_S = self.get_local_action(x, y, z, t)
                        
                        proposal = old_val + np.random.uniform(-self.delta, self.delta)
                        self.phi[x, y, z, t] = proposal
                        new_S = self.get_local_action(x, y, z, t)
                        
                        dS = new_S - old_S
                        if dS > 0 and np.random.rand() > np.exp(-dS):
                            self.phi[x, y, z, t] = old_val
                        else:
                            accepted += 1
        return accepted / (self.L**4)

    def calculate_correlation(self, max_r):
        """距離 r に対する相関関数 G(r) を計算 (軸方向に平均)"""
        correlations = np.zeros(max_r)
        for r in range(max_r):
            sum_corr = 0
            for axis in range(4):
                rolled = np.roll(self.phi, -r, axis=axis)
                sum_corr += np.mean(self.phi * rolled)
            correlations[r] = sum_corr / 4.0
        return correlations
