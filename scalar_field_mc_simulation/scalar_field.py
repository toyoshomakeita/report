import numpy as np

class ScalarField2D:
    def __init__(self, L, m, delta=0.5):
        """
        2次元実スカラー場の格子シミュレーション・クラス
        
        Args:
            L (int): 格子サイズ (L x L)
            m (float): 質量
            delta (float): メトロポリス更新のステップサイズ
        """
        self.L = L
        self.m = m
        self.delta = delta
        self.phi = np.zeros((L, L)) # 冷たいスタート (All zeros)
        
    def get_local_action(self, i, j):
        """指定された点 (i, j) に依存する局所的な作用を計算する"""
        # 周期境界条件
        phi_curr = self.phi[i, j]
        phi_neighbors = (
            self.phi[(i + 1) % self.L, j] +
            self.phi[(i - 1) % self.L, j] +
            self.phi[i, (j + 1) % self.L] +
            self.phi[i, (j - 1) % self.L]
        )
        
        # 運動項: 1/2 * sum_{mu} (phi(x+mu) - phi(x))^2
        # 展開すると: sum_{mu} (phi(x)^2 - phi(x)phi(x+mu) + 1/2 phi(x+mu)^2)
        # phi(x) を含む項のみを抽出: 2 * phi(x)^2 - phi(x) * sum_{neighbor} phi(neighbor)
        kinetic = 2.0 * phi_curr**2 - phi_curr * phi_neighbors
        
        # 質量項: 1/2 * m^2 * phi(x)^2
        mass_term = 0.5 * (self.m**2) * (phi_curr**2)
        
        return kinetic + mass_term

    def metropolis_sweep(self):
        """全格子点に対して1回ずつメトロポリス更新を試みる (1 Sweep)"""
        accepted = 0
        for i in range(self.L):
            for j in range(self.L):
                old_phi = self.phi[i, j]
                old_action = self.get_local_action(i, j)
                
                # 新しい値を提案
                proposal = old_phi + np.random.uniform(-self.delta, self.delta)
                self.phi[i, j] = proposal
                new_action = self.get_local_action(i, j)
                
                dS = new_action - old_action
                
                # メトロポリス判定
                if dS > 0 and np.random.rand() > np.exp(-dS):
                    # 棄却
                    self.phi[i, j] = old_phi
                else:
                    accepted += 1
        
        return accepted / (self.L**2)

    def calculate_total_action(self):
        """格子全体の全作用を計算する"""
        action = 0.0
        for i in range(self.L):
            for j in range(self.L):
                # 運動項 (一方向のみ計算して重複を防ぐ)
                phi_curr = self.phi[i, j]
                diff_x = self.phi[(i + 1) % self.L, j] - phi_curr
                diff_y = self.phi[i, (j + 1) % self.L] - phi_curr
                
                action += 0.5 * (diff_x**2 + diff_y**2)
                action += 0.5 * (self.m**2) * (phi_curr**2)
        return action

    def measure_phi2(self):
        """phi^2 の平均値を計算する"""
        return np.mean(self.phi**2)

    def calculate_correlation_function(self, max_r):
        """
        距離 r に対する相関関数 G(r) = <phi(0)phi(r)> を計算する
        簡単のため、x方向の相関のみを計算して平均をとる
        """
        correlations = np.zeros(max_r)
        for r in range(max_r):
            # phi(x, y) と phi(x+r, y) の積の平均
            rolled_phi = np.roll(self.phi, -r, axis=0)
            correlations[r] = np.mean(self.phi * rolled_phi)
        return correlations
