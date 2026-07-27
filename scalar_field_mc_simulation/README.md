# 2次元実スカラー場の格子モンテカルロシミュレーション

## 概要

2次元ユークリッド空間上の実スカラー場 \(\phi(x)\) を格子場理論の手法でシミュレーションし、
\(\langle \phi^2 \rangle\) および相関関数 \(G(r) = \langle \phi(0)\phi(r) \rangle\) を測定して
質量 \(m\) を推定する。

## 物理的背景

### ユークリッド作用

2次元連続空間におけるスカラー場の作用：

\[ S[\phi] = \int d^2x \left[ \frac{1}{2} \sum\_{\mu=1}^2 (\partial\_\mu \phi)^2 + \frac{1}{2} m^2 \phi^2 \right] \]

### 格子上の局所作用

離散化・局所化すると（隣接点は上下左右の4点）：

\[ S\_{\text{local}}(\phi\_n) = \left( 2 + \frac{1}{2} m^2 \right) \phi\_n^2 - \phi\_n \sum\_{\text{neighbors}} \phi\_{\text{neighbor}} \]

### 相関関数

経路積分による期待値は統計力学的なアンサンブル平均と等価であり、メトロポリス法で
重み \(e^{-S}\) に従うサンプルを生成することで：

\[ G(r) = \langle \phi(0)\phi(r) \rangle \approx \frac{1}{N} \sum\_{i=1}^N \phi^{(i)}(0)\phi^{(i)}(r) \]

遠方で \(G(r) \sim e^{-mr}\) の指数減衰を示すため、\(\ln G(r)\) の傾きから有効質量を推定する。

## 実装

### ファイル構成

| ファイル | 説明 |
| --- | --- |
| `scalar_field.py` | `ScalarField2D` クラス：2次元格子場の定義、メトロポリス更新、全作用・\(\phi^2\)・相関関数の測定 |
| `main.py` | シミュレーション実行スクリプト |
| `correlation_results.png` | 相関関数のプロット図（3パネル：場の配置・\(G(r)\)・対数プロット） |
| `simulation_results.png` | シミュレーション結果の図 |

### コード構造

* **`ScalarField2D(L, m, delta)`**: 2次元配列 `phi[L][L]` を保持
  + `get_local_action(i, j)`: 指定点の局所作用（周期境界条件）
  + `metropolis_sweep()`: 全 \(L^2\) 点を1スイープ、採択率を返す
  + `calculate_total_action()`: 全作用を計算（熱化の経過確認用）
  + `measure_phi2()`: \(\langle \phi^2 \rangle\) を測定
  + `calculate_correlation_function(max_r)`: x方向の相関関数を計算
* **`main.py`**: \(L=20\)、\(m=0.5\)、熱化500sweep、測定1000sweep（interval=10）、\(r=1\) から \(L/4\) の範囲で線形フィット

## 結果と考察

* **入力質量**: \(m = 0.5\)
* **測定量**: \(\langle \phi^2 \rangle\) および相関関数 \(G(r)\)
* **質量推定**: \(\ln G(r)\) の線形フィットから有効質量 \(m\_{\text{eff}}\) を算出

**注：**このプロジェクトは2次元（\(L=20\)、\(20^2=400\)点）であり、
`scalar_field_mc_4d` の4次元シミュレーションとは次元数・格子サイズ・物理的セットアップが異なる。
2次元では相関関数は遠方で \(G(r) \sim \frac{1}{\sqrt{r}} e^{-mr}\) の振る舞いを示す（4次元の \(1/r^2\) とは異なる）。

## 参考文献・HTMLドキュメント

| ファイル | 内容 |
| --- | --- |
| `explanation.html` | 理論解説（作用・離散化・局所作用・モンテカルロ法） |
| `explanation.md` | 同内容のMarkdown版 |
| `correlation_explanation.html` | 相関関数と経路積分の詳細解説（重点サンプリング・メトロポリス法） |

2026年6月作成 | 2次元格子場理論シミュレーション