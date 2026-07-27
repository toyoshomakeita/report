# 4次元スカラー場：作用の数学的導出プロセス

このノートでは、シミュレーションコードで使用している「局所作用」の式が、物理学の基本的な定義からどのように導かれるかを、1ステップずつ式変形を追って解説します。

## 1. 出発点：ミンコフスキー空間の作用

まず、私たちが住む4次元空間（時間1次元＋空間3次元）でのスカラー場の作用 \(S\_M\) は次のように定義されます：

\[ S\_M = \int d^4x \left[ \frac{1}{2} \eta^{\mu\nu} \partial\_\mu \phi \partial\_\nu \phi - \frac{1}{2} m^2 \phi^2 \right] \]

ここで、\(\eta^{\mu\nu} = \text{diag}(1, -1, -1, -1)\) です。書き下すと：

\[ S\_M = \int dt d^3x \left[ \frac{1}{2} (\partial\_t \phi)^2 - \frac{1}{2} (\nabla \phi)^2 - \frac{1}{2} m^2 \phi^2 \right] \]

ここでの問題：経路積分 \(\int \mathcal{D}\phi e^{iS\_M}\) の指数部が「虚数 \(i\)」であるため、確率として扱えません。

## 2. ステップ1：ウィック回転 (Wick Rotation)

時間を虚数化します：\(t \to -i\tau\)。このとき、微分の連鎖律より：

\[ \partial\_t = \frac{\partial \tau}{\partial t} \partial\_\tau = i \partial\_\tau \]

また、積分測度は \(dt \to -id\tau\) となります。これらを \(S\_M\) に代入します：

\[ i S\_M = i (-i) \int d\tau d^3x \left[ \frac{1}{2} (i \partial\_\tau \phi)^2 - \frac{1}{2} (\nabla \phi)^2 - \frac{1}{2} m^2 \phi^2 \right] \]
\[ i S\_M = \int d^4x\_E \left[ -\frac{1}{2} (\partial\_\tau \phi)^2 - \frac{1}{2} (\nabla \phi)^2 - \frac{1}{2} m^2 \phi^2 \right] \]

ここで、\(e^{i S\_M} = e^{-S\_E}\) と定義すると、ユークリッド作用 \(S\_E\) が得られます：

\[ S\_E = \int d^4x\_E \left[ \frac{1}{2} \sum\_{\mu=1}^4 (\partial\_\mu \phi)^2 + \frac{1}{2} m^2 \phi^2 \right] \]

全てがプラスの項になり、重みが実数（確率）になりました！

## 3. ステップ2：格子への離散化 (Discretization)

連続な積分を、間隔 \(a=1\) の格子上の和に置き換えます。微分項に注目しましょう：

\[ \int d^4x \frac{1}{2} \sum\_{\mu=1}^4 (\partial\_\mu \phi)^2 \quad \longrightarrow \quad \sum\_{x} \sum\_{\mu=1}^4 \frac{1}{2} (\phi(x+\hat{\mu}) - \phi(x))^2 \]

この二乗を展開します：

\[ \frac{1}{2} (\phi(x+\hat{\mu})^2 - 2\phi(x)\phi(x+\hat{\mu}) + \phi(x)^2) \]

格子全体で和をとると、全ての点 \(x\) について \(\phi(x+\hat{\mu})^2\) と \(\phi(x)^2\) が現れるため、結局同じ項が2倍ずつカウントされることになります。

## 4. ステップ3：局所作用の抽出 (Local Action)

メトロポリス法では、特定の1点 \(x\_0\) の値だけを更新します。したがって、**「\(\phi(x\_0)\) を含んでいる項」**だけを取り出せば十分です。

### 1. 運動項からの寄与

\(\sum\_{x, \mu} (\phi(x+\hat{\mu}) - \phi(x))^2\) の中で、\(x\_0\) が関わるのは以下の2つのケースです：

* \(x = x\_0\) のとき： \((\phi(x\_0+\hat{\mu}) - \phi(x\_0))^2\)
* \(x+\hat{\mu} = x\_0\) のとき： \((\phi(x\_0) - \phi(x\_0-\hat{\mu}))^2\)

これらを展開して \(\phi(x\_0)\) に関する項（2乗の項と隣接点との積の項）を抽出すると：

\[ \text{from } x=x\_0: \quad \frac{1}{2} \sum\_{\mu=1}^4 (\phi(x\_0)^2 - 2\phi(x\_0)\phi(x\_0+\hat{\mu})) \]
\[ \text{from } x+\hat{\mu}=x\_0: \quad \frac{1}{2} \sum\_{\mu=1}^4 (\phi(x\_0)^2 - 2\phi(x\_0)\phi(x\_0-\hat{\mu})) \]

これらを足し合わせます：

\[ \text{Kinetic Part} = \frac{1}{2} \left[ 4\phi(x\_0)^2 + 4\phi(x\_0)^2 - 2\phi(x\_0) \sum\_{\mu=1}^4 (\phi(x\_0+\hat{\mu}) + \phi(x\_0-\hat{\mu})) \right] \]
\[ = 4\phi(x\_0)^2 - \phi(x\_0) \sum\_{\text{neighbor}} \phi(x\_{\text{neighbor}}) \]

### 2. 質量項からの寄与

これは単純です：

\[ \text{Mass Part} = \frac{1}{2} m^2 \phi(x\_0)^2 \]

## 5. 結論：完成した局所作用の式

以上の項を全て合算すると、コード内で実装されている式が完成します：

\[ S\_{\text{local}}(x\_0) = 4 \phi(x\_0)^2 - \phi(x\_0) \sum\_{n=1}^8 \phi(x\_n) + \frac{1}{2} m^2 \phi(x\_0)^2 \]

コード内の `scalar_field_4d.py` の `get_local_action` 関数：

```
kinetic = 4.0 * phi_curr**2 - phi_curr * phi_neighbors
mass_term = 0.5 * (self.m**2) * (phi_curr**2)
return kinetic + mass_term
```

**まとめ：** 係数の「4」は次元数 \(d=4\) から来ており、隣接点の和が「8点」なのは 4方向×2（前後）だからです。