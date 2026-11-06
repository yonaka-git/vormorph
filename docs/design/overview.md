# 設計概要：支配方程式と離散化

この文書は、vormorph の支配方程式、ボロノイ C 格子上の離散式、1ステップの計算手順をまとめたものです。設計判断の理由は [ADR](adr/README.md) に記録しています。

> **状態**：流れの半陰解法（4.1〜4.3節）、wet/dry、移流（4.4節）、Perot 法による再構成（4.5節）、境界条件（6章。直線水路の周期境界を含む）、時間刻みの決め方（7章）、掃流砂と河床の更新（4.6節、4.8節）は実装済みです。浮遊砂（4.7節と、4.8節の浮遊砂の項）はまだ実装していません。実装の過程で式や手順を変えた場合は、この文書と ADR を同時に更新します。末尾の「[未決事項](#8-未決事項)」も参照してください。

## 目次

1. [対象と前提](#1-対象と前提)
2. [支配方程式](#2-支配方程式)
3. [格子と変数配置](#3-格子と変数配置)
4. [離散式](#4-離散式)
5. [1ステップの手順](#5-1ステップの手順)
6. [境界条件](#6-境界条件)
7. [時間刻みの制約](#7-時間刻みの制約)
8. [未決事項](#8-未決事項)

## 1. 対象と前提

- 対象は、堰や床止めを含む河道の平面二次元の流れと河床変動です。跳水や射流が局所的に生じます。
- 流れは水深平均の浅水流方程式で表し、静水圧を仮定します。
- 流砂は均一粒径の掃流砂と浮遊砂として扱います。混合粒径（平野の交換層モデル）は将来の拡張です。
- 単位はすべて SI 単位とし、浮動小数点数は `f64` を使います。

## 2. 支配方程式

### 2.1 流れ

水位を $\eta$、河床高を $z_b$、水深を $h = \eta - z_b$、水深平均流速を $\mathbf{u}$、単位幅流量を $\mathbf{q} = h\mathbf{u}$ とします。

連続式：

$$
\frac{\partial \eta}{\partial t} + \nabla\cdot(h\mathbf{u}) = 0
$$

運動量式（摩擦はマニング式）：

$$
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u} = -g\nabla\eta - \frac{g n^2 \lvert\mathbf{u}\rvert}{h^{4/3}}\mathbf{u}
$$

水平渦粘性は当面省略します。移流項は、跳水で運動量が保存されるように、次の保存形に書き直してから離散化します（Stelling & Duinmeijer, 2003）。

$$
(\mathbf{u}\cdot\nabla)\mathbf{u} = \frac{1}{h}\left[\nabla\cdot(\mathbf{q}\otimes\mathbf{u}) - \mathbf{u}\,\nabla\cdot\mathbf{q}\right]
$$

### 2.2 掃流砂

水中比重を $s = \rho_s/\rho - 1$、粒径を $d$ とし、無次元掃流力と無次元掃流砂量を次のように定義します。

$$
\tau_\ast = \frac{u_\ast^2}{s g d}, \qquad u_\ast^2 = \frac{g n^2 \lvert\mathbf{u}\rvert^2}{h^{1/3}}, \qquad q_b = \sqrt{s g d^3}\;\Phi(\tau_\ast)
$$

$\Phi$ は Meyer-Peter & Müller (1948) の式と芦田・道上 (1972) の式から、設定で選びます（$\tau_\ast \le \tau_{\ast c}$ では 0）。

$$
\Phi_{\mathrm{MPM}} = 8\left(\tau_\ast - \tau_{\ast c}\right)^{3/2}, \qquad \Phi_{\mathrm{AM}} = 17\,\tau_{\ast e}^{3/2}\left(1 - \frac{\tau_{\ast c}}{\tau_\ast}\right)\left(1 - \sqrt{\frac{\tau_{\ast c}}{\tau_\ast}}\right)
$$

$\tau_{\ast e}$ は有効掃流力です。当面は $\tau_{\ast e} = \tau_\ast$ とします。

掃流砂の向きは流速の向きを基本とし、河床勾配による補正を長谷川の式（長谷川, 1984。式の形は長谷川, 1983 と Nays2DH マニュアルで確認）で加えます。任意の方向 $\mathbf{n}$ の成分は次のとおりです。

$$
\mathbf{q}_b\cdot\mathbf{n} = q_b\left(\frac{\mathbf{u}_b\cdot\mathbf{n}}{\lvert\mathbf{u}_b\rvert} - \gamma\,\frac{\partial z_b}{\partial n}\right), \qquad \gamma = \sqrt{\frac{\tau_{\ast c}}{\mu_s \mu_k \tau_\ast}}
$$

$\mu_s$ は静止摩擦係数、$\mu_k$ は動摩擦係数です。$\mathbf{u}_b$ は河床近傍の流速で、二次流を考えない当面の実装では水深平均流速 $\mathbf{u}$ で代用します。

### 2.3 浮遊砂

水深平均濃度を $C$ とし、$hC$ を保存量とする移流拡散式を解きます。

$$
\frac{\partial (hC)}{\partial t} + \nabla\cdot(\mathbf{q}\,C) = \nabla\cdot(h\varepsilon\nabla C) + E - w_f c_b
$$

- $E$ は浮上量で、板倉・岸の式（Itakura & Kishi, 1980）または Garcia & Parker (1991) の式を使います。
- $w_f$ は沈降速度で、Rubey (1933) の式から求めます。
- $c_b$ は河床近傍の濃度で、$c_b = r\,C$ と表します。たとえば鉛直分布に指数分布を仮定すると $r = \beta/(1 - e^{-\beta})$、$\beta = 6 w_f/(\kappa u_\ast)$ です。$r$ の具体形は浮上量式に合わせて実装時に決めます。
- 水平拡散係数は $\varepsilon = \alpha\, u_\ast h$ とします。

### 2.4 河床（Exner 式）

空隙率を $\lambda$ とします。

$$
(1-\lambda)\frac{\partial z_b}{\partial t} + \nabla\cdot\mathbf{q}_b + E - w_f c_b = 0
$$

地形変化の時間スケールを縮めるため、河床の更新量に MORFAC $M$ を掛けます（流れと浮遊砂の式には掛けません）。

## 3. 格子と変数配置

![ボロノイ C 格子の変数配置](../images/c-grid.svg)

ボロノイ格子では、隣り合う母点を結ぶ線分と、両者を隔てる辺が直交します。この性質を使い、変数を C 格子（スタガード）に配置します。

- セル中心（母点）：水位 $\eta_i$、河床高 $z_{b,i}$、浮遊砂 $(hC)_i$
- 辺：法線流速 $u_j$

| 記号 | 意味 | コード上の名前 |
|---|---|---|
| $i$, $j$ | セル番号、辺番号 | |
| $1(j)$, $2(j)$ | 辺 $j$ を挟む2つのセル | `cells_on_edge[j][0]`, `[1]` |
| $k(i,j)$ | 辺 $j$ を挟んでセル $i$ の隣にあるセル | `cells_on_cell` |
| $\partial i$ | セル $i$ を囲む辺の集合 | `edges_on_cell` |
| $s_{i,j}$ | 辺 $j$ の法線がセル $i$ から外向きなら $+1$、内向きなら $-1$ | `edge_sign_on_cell` |
| $A_i$ | セル面積 | `area_cell` |
| $L_j$ | 辺の長さ | `dv_edge` |
| $\delta_j$ | 母点間の距離 | `dc_edge` |
| $\mathbf{n}_j$ | 辺の単位法線（$1(j)$ から $2(j)$ へ向く） | `normal_edge` |
| $\mathbf{x}_i$, $\mathbf{x}_j$ | 母点の座標、辺の中点の座標 | `x_cell`, `x_edge` など |
| $H_j$ | 辺の水深 | |
| $\theta$ | θ 法の重み（0.55〜0.6） | |

符号の規約は次のとおりです。

- $u_j > 0$ は $1(j)$ から $2(j)$ へ向かう流れを表します。
- セル $i$ から外へ出る流量は $\sum_{j\in\partial i} s_{i,j} L_j H_j u_j$ です。
- 境界辺では $2(j)$ を領域外（`NONE = usize::MAX`）とし、法線は領域の外を向きます。
- 周期境界（6章）の継ぎ目の辺は内部の辺です。$2(j)$ の母点は、$1(j)$ の側から見た位置（`x_cell` に `edge_shift` を足した位置）で扱います。辺の中点と頂点は $1(j)$ の側の座標で持ちます。

MPAS-Ocean の `edgeSignOnCell` はこれと逆符号（`cellsOnEdge(1)` 側を $-1$）なので、MPAS 形式のデータを読み込むときは符号を反転します。

## 4. 離散式

### 4.1 運動量式（辺）

移流を陽的に、水位勾配を θ 法で、摩擦を陰的に扱います（Casulli & Walters, 2000）。

$$
u_j^{n+1} = F u_j^n - \frac{g\Delta t}{\delta_j}\left[\theta\left(\eta_{2(j)}^{n+1} - \eta_{1(j)}^{n+1}\right) + (1-\theta)\left(\eta_{2(j)}^{n} - \eta_{1(j)}^{n}\right)\right] - \Delta t\,\gamma_j\,u_j^{n+1}
$$

$F u_j^n$ は移流を評価した後の流速（4.4節）です。$\gamma_j$ は摩擦係数で、時刻 $n$ の値で線形化します。

$$
\gamma_j = \frac{g n^2 \lvert\mathbf{u}\rvert_j^n}{\left(H_j^n\right)^{4/3}}
$$

$u_j^{n+1}$ について解くと、次の形になります。

$$
u_j^{n+1} = G_j - \frac{\theta g \Delta t}{\delta_j\left(1 + \Delta t\,\gamma_j\right)}\left(\eta_{2(j)}^{n+1} - \eta_{1(j)}^{n+1}\right)
$$

$$
G_j = \frac{1}{1 + \Delta t\,\gamma_j}\left[F u_j^n - (1-\theta)\frac{g\Delta t}{\delta_j}\left(\eta_{2(j)}^{n} - \eta_{1(j)}^{n}\right) + \Delta t\, g I_p\, n_{x,j}\right]
$$

$g I_p n_{x,j}$ は周期境界（6章）の体積力です。$I_p$ は、河床高と水位から除いた縦断勾配で、周期境界でなければ 0 です。勾配を除いた水位の差は $\partial\eta/\partial n = (\eta_2 - \eta_1)/\delta_j - I_p n_{x,j}$ なので、θ 法の両方の時刻から同じ定数が出ます（[ADR-0014](adr/0014-periodic-channel.md)）。

辺の水深 $H_j^n$ は、風上側のセルの水位から辺の河床高を引いた値とします（Stelling & Duinmeijer, 2003）。

$$
H_j^n = \max\left(\eta_{\mathrm{up}(j)}^n - z_{b,j},\ 0\right), \qquad z_{b,j} = \max\left(z_{b,1(j)},\ z_{b,2(j)}\right)
$$

辺の河床高に両側の高い方を使うと、床止めのような段差で、天端を越える流れの水深を正しく表せます。細部は次のとおりです（[ADR-0006](adr/0006-wet-dry-and-newton.md)）。

- 風上側は $u_j^n$ の向きで決めます。$u_j^n = 0$ のときは、両側の水位の高い方を使います。
- $H_j^n$ が閾値 `dry_depth` 未満の辺は乾いているとみなし、$H_j^n = 0$、$u_j^{n+1} = 0$ とします（水位の連立方程式からも外れます）。セルの体積関数（4.2節）には閾値を入れません。
- 周期境界（6章）では、河床高と水位を縦断勾配 $I_p$ を除いた値で持つので、セル 2 の値から除いた勾配の分 $I_p(x_2 - x_1)$ を引いて、実際の高さで風上と高い方の河床を選びます（[ADR-0015](adr/0015-periodic-edge-depth.md)）。勾配を除いたまま選ぶと、河床が下流へ上がる所で下流側のセルが選ばれ（風下）、河床に 2Δx の振動が育ちます。
- 静水（$\eta$ が一定、$u = 0$）では、水面より高い乾いたセル（水位を河床高に揃えている。4.3節）と接する辺の水深は $z_b - z_b = 0$ になります。したがって、wet/dry の境界で偽の流れは生じません（検証ケース V03）。

摩擦係数の $\lvert\mathbf{u}\rvert_j^n$ は、法線成分 $u_j^n$ と、両側のセルの Perot 流速（4.5節）の平均の接線成分 $\bar{\mathbf{u}}_j\cdot\mathbf{t}_j$ から求めます。

$$
\lvert\mathbf{u}\rvert_j^n = \sqrt{\left(u_j^n\right)^2 + \left(\bar{\mathbf{u}}_j\cdot\mathbf{t}_j\right)^2}, \qquad \bar{\mathbf{u}}_j = \tfrac{1}{2}\left(\mathbf{u}_{1(j)}^n + \mathbf{u}_{2(j)}^n\right)
$$

Perot 法は一様流で厳密なので（V02）、一様流ではこれも元の速さに厳密に一致します。境界辺では内側のセルの流速を使います。

### 4.2 連続式（セル）

セル $i$ の水の体積を $V_i(\eta_i)$ とします。セル内の河床高が一定なら、次のようになります。

$$
V_i(\eta_i) = A_i \max\left(\eta_i - z_{b,i},\ 0\right)
$$

連続式は、辺の流量に θ 平均の流速を使って次のように離散化します。

$$
V_i\left(\eta_i^{n+1}\right) = V_i\left(\eta_i^{n}\right) - \Delta t\sum_{j\in\partial i} s_{i,j}\, L_j\, H_j^n\left[\theta\, u_j^{n+1} + (1-\theta)\, u_j^n\right] + \Delta t\, Q_i
$$

$Q_i$ は上流端の流量境界からの流入です。この離散式はセルごとに体積を厳密に保存し、水位方程式の解の精度（Newton 法と PCG の収束判定値）の範囲で全体の質量も保存されます。

### 4.3 水位の連立方程式

4.1節の $u_j^{n+1}$ を4.2節に代入します。$s_{i,j}\left(\eta_{2(j)} - \eta_{1(j)}\right) = \eta_{k(i,j)} - \eta_i$ を使うと、水位について次の連立方程式が得られます。

$$
V_i\left(\eta_i^{n+1}\right) + \sum_{j\in\partial i} w_j\left(\eta_i^{n+1} - \eta_{k(i,j)}^{n+1}\right) = b_i, \qquad w_j = \frac{g\,\theta^2 \Delta t^2\, L_j\, H_j^n}{\delta_j\left(1 + \Delta t\,\gamma_j\right)}
$$

$$
b_i = V_i\left(\eta_i^{n}\right) - \Delta t\sum_{j\in\partial i} s_{i,j}\, L_j\, H_j^n\left[\theta\, G_j + (1-\theta)\, u_j^n\right] + \Delta t\, Q_i
$$

ベクトルで書くと $\mathbf{V}(\boldsymbol{\eta}) + \mathsf{T}\boldsymbol{\eta} = \mathbf{b}$ です。$\mathsf{T}$ は重み $w_j \ge 0$ を持つグラフラプラシアン型の対称行列で、半正定値です。

$\mathbf{V}$ は wet/dry の境目（$\eta_i = z_{b,i}$）で折れる区分線形の関数なので、Newton 法で解きます（Casulli, 2009; Casulli & Zanolli, 2012）。

$$
\left[\mathsf{P}^{(m)} + \mathsf{T}\right]\boldsymbol{\eta}^{(m+1)} = \mathbf{b} - \mathbf{V}\left(\boldsymbol{\eta}^{(m)}\right) + \mathsf{P}^{(m)}\boldsymbol{\eta}^{(m)}
$$

- $\mathsf{P}^{(m)}$ は $\mathbf{V}$ の微分で、$\eta_i^{(m)} \ge z_{b,i}$ のセルでは $A_i$、それ以外のセルでは 0 の対角行列です。折れ点（$\eta_i = z_{b,i}$）では右側の微分 $A_i$ をとります。
- $\mathsf{P} + \mathsf{T}$ は、各連結成分に $\mathsf{P}_i > 0$ のセルか水位境界が1つでもあれば対称正定値です。各段の線形系は、行列を組み立てずに前処理付き共役勾配法（Jacobi 前処理）で解きます。
- 行列とベクトルの積はセル単位のループで計算し、各セルが自分の成分だけを書き込むようにします。将来 rayon で並列化するときに、書き込みの競合が起きないようにするためです。

実装では、同じ式を修正量 $\boldsymbol{\delta}^{(m)} = \boldsymbol{\eta}^{(m+1)} - \boldsymbol{\eta}^{(m)}$ の形に書き直して解きます（[ADR-0006](adr/0006-wet-dry-and-newton.md)）。

$$
\left[\mathsf{P}^{(m)} + \mathsf{T}\right]\boldsymbol{\delta}^{(m)} = \mathbf{r}^{(m)}, \qquad \mathbf{r}^{(m)} = \mathbf{b} - \mathbf{V}\left(\boldsymbol{\eta}^{(m)}\right) - \mathsf{T}\boldsymbol{\eta}^{(m)}
$$

- PCG の相対判定値は $\mathbf{r}^{(m)}$ に対してかかります。Newton 反復が進むほど解の精度が上がり、最終的な残差（水の収支の誤差になる）は収束判定値よりずっと小さくなります。
- 収束判定は $\max_i \lvert r_i\rvert / A_i \le$ `newton_tol` です。乾いたセルも含めて判定するので、負の体積が必要な解は、収束の失敗として検出されます。

**初期値と収束**：初期値は $\boldsymbol{\eta}^{(0)} = \boldsymbol{\eta}^n$ とします。各ステップの最後に乾いたセルの水位を河床高に揃える（$\eta_i \leftarrow \max(\eta_i, z_{b,i})$。$V_i$ は変わりません）ので、全セルで $\eta_i^n \ge z_{b,i}$ です。したがって、最初の段の係数行列は $\mathrm{diag}(A) + \mathsf{T}$ で、対称正定値です。

収束は次のように示せます。$\mathbf{f}(\boldsymbol{\eta}) = \mathbf{V}(\boldsymbol{\eta}) + \mathsf{T}\boldsymbol{\eta} - \mathbf{b}$、解を $\boldsymbol{\eta}^\ast$、$\mathsf{J}^{(m)} = \mathsf{P}^{(m)} + \mathsf{T}$ とします。

1. $V_i$ は凸なので、任意の $\boldsymbol{\eta}$ について $\mathbf{f}(\boldsymbol{\eta}) \ge \mathbf{f}(\boldsymbol{\eta}^{(m)}) + \mathsf{J}^{(m)}(\boldsymbol{\eta} - \boldsymbol{\eta}^{(m)})$（成分ごと）です。
2. $\boldsymbol{\eta} = \boldsymbol{\eta}^\ast$ とおくと $\mathsf{J}^{(m)}(\boldsymbol{\eta}^{(m+1)} - \boldsymbol{\eta}^\ast) \ge \mathbf{0}$ が得られます。$\mathsf{J}^{(m)}$ は非対角成分が 0 以下の正定値行列（M 行列）なので逆行列の成分は非負で、$\boldsymbol{\eta}^{(m+1)} \ge \boldsymbol{\eta}^\ast$ です。
3. $\boldsymbol{\eta} = \boldsymbol{\eta}^{(m+1)}$ とおくと $\mathbf{f}(\boldsymbol{\eta}^{(m+1)}) \ge \mathbf{0}$ なので、次の修正量は 0 以下で、2段目以降の反復値は単調に減少します。
4. 反復値は解以上にあるので、解で wet のセルは wet のままで、係数行列は正定値を保ちます。$\mathbf{V}$ は区分線形なので、wet/dry の判定が変わらなくなった段で収束します。

これは Casulli (2009)、Casulli & Zanolli (2012) と同じ議論です。壁で囲まれ、乾いた辺（$w_j = 0$）で区切られた領域では、辺の流量が打ち消し合うので $\sum_i b_i = \sum_i V_i^n \ge 0$ です。したがって、水があれば解でも wet のセルが残り、解が存在します。

$\boldsymbol{\eta}^{n+1}$ が求まったら、4.1節の式で $u_j^{n+1}$ を更新します。

### 4.4 移流（保存型風上差分）

移流は陽的な保存型風上差分で評価します（Stelling & Duinmeijer, 2003; Kramer & Stelling, 2008）。決め方の理由は [ADR-0008](adr/0008-advection-discretization.md) と [ADR-0009](adr/0009-advection-face-velocity.md) にあります。

2.1節の保存形 $h^{-1}\left[\nabla\cdot(\mathbf{q}\otimes\mathbf{u}) - \mathbf{u}\,\nabla\cdot\mathbf{q}\right]$ をセル $i$ で評価します。運動量の流束の風上の流速は、流出側の辺では $\mathbf{u}_i$ です。そこで、連続式の項に辺の流速 $u_j$ を使うと、流入側の辺だけが残ります。

まず、セルごとに、流入量 $S_i$ と、流れ込む運動量 $\mathbf{P}_i$ を求めます。

$$
S_i = \sum_{l\in\partial i,\ s_{i,l} q_l < 0} \lvert L_l q_l\rvert, \qquad \mathbf{P}_i = \frac{1}{A_i}\sum_{l\in\partial i,\ s_{i,l} q_l < 0} \lvert L_l q_l\rvert\,\mathbf{u}_{k(i,l)}, \qquad q_l = H_l^n u_l^n
$$

$\mathbf{u}_k$ は4.5節の Perot 法で再構成したセル中心の流速です。境界辺から流れ込む水の流速は、その辺の法線流速 $u_l\mathbf{n}_l$ とします（流量境界の流入）。

次に、時刻 $n$ の同じ流量で連続式を陽的に進めたセルの水深 $h^\ast_i$ を求めます。辺の流速を、両側のセルに流れ込む水の流速（法線成分）$\bar u_j$ へ、局所クーラン数 $\sigma_j$ の割合で近づけます。

$$
h^\ast_i = h_i^n - \frac{\Delta t}{A_i}\sum_{l\in\partial i} s_{i,l}\, L_l\, q_l, \qquad \sigma_j = \Delta t\,\frac{S_{1(j)}/A_{1(j)} + S_{2(j)}/A_{2(j)}}{h^\ast_{1(j)} + h^\ast_{2(j)}}, \qquad \bar u_j = \frac{\left(\mathbf{P}_{1(j)} + \mathbf{P}_{2(j)}\right)\cdot\mathbf{n}_j}{S_{1(j)}/A_{1(j)} + S_{2(j)}/A_{2(j)}}
$$

$$
F u_j^n = u_j^n + \min(\sigma_j, 1)\left(\bar u_j - u_j^n\right)
$$

- **重みは両側とも 1/2 です。** ボロノイ格子の辺は、母点を結ぶ線分の垂直二等分線だからです。
- **跳水で Bélanger の式が成り立ちます。** 定常流では $h^\ast = h$ です。水位勾配の項 $g\,\bar h_j\,\partial\eta/\partial n$（$\bar h_j = (h_1 + h_2)/2$）と合わせると、1次元の流れでは、各辺の式の和で $q u + g h^2/2$ が打ち消し合います（連続式の項に $u_j$ を使うことで加わる項は流速の2階差分で、和では両端の差だけが残ります）。そのため、跳水の両側で Rankine–Hugoniot の関係（Bélanger の式）が離散式の上で成り立ちます。
- **$h^\ast$ で割ると、陽的な有限体積法の更新と一致します。** セルの運動量 $V\mathbf{u}$ と体積 $V$ を同じ流量で陽的に進めた形です。7章の時間刻みの条件が成り立てば $\sigma_j \le 1$ で、$F u_j^n$ は $u_j^n$ と $\bar u_j$ の凸結合になります。濡れたばかりの浅いセルや段波の先端でも発散しません。
- **格子規模の辺の流速の成分も運び去られます。** 辺の流速そのものを上流の流速へ近づけるためです。六角形配置では、辺の数がセル数の約3倍あり、セル中心の流速に現れない格子規模の辺の流速の成分が生じます。セル中心の流速だけで増分を作ると、この成分が跳水などにたまります（ADR-0009）。
- **安全策**：$\sigma_j > 1$ の辺（時間刻みの条件を破ったとき）では、重みを 1 にします（$F u_j^n = \bar u_j$）。
- **水位境界の辺**は、内側のセルの値だけで評価します。壁と流量境界の辺には運動量式がないので、移流も評価しません。
- 1次精度の風上差分です。2次精度化（制限関数付き）と、加速流（堰の天端へ向かう流れなど）でのエネルギー水頭保存形への切り替えは、未決事項（8章）です。
- 設定 `[hydro] advection = "none"` で、移流項を省けます（$F u_j^n = u_j^n$。検証で移流の効果を比べるため）。

### 4.5 流速の再構成（Perot 法）

セル中心の流速ベクトルを、辺の法線流速から次の式で再構成します（Perot, 2000）。

$$
\mathbf{u}_i = \frac{1}{A_i}\sum_{j\in\partial i} s_{i,j}\, L_j\, u_j\left(\mathbf{x}_j - \mathbf{x}_i\right)
$$

ここで $\mathbf{x}_j$ は辺の中点です。一様流 $\mathbf{u}$ に対しては、発散定理から次の式が成り立ちます。

$$
\oint_{\partial\Omega_i} (\mathbf{u}\cdot\mathbf{n})\left(\mathbf{x} - \mathbf{x}_i\right) ds = \int_{\Omega_i} \nabla\cdot\left[\mathbf{u}\otimes\left(\mathbf{x} - \mathbf{x}_i\right)\right] dA = A_i\,\mathbf{u}
$$

辺の上では被積分関数が $\mathbf{x}$ の1次式です。したがって辺の中点で評価した上の和はこの線積分に一致し、再構成は一様流で厳密になります（検証ケース V02）。周期境界の継ぎ目の辺では、$\mathbf{x}_j$ をセルの側の座標で使います。

### 4.6 掃流砂（辺の法線成分）

掃流砂量は、辺 $j$ の法線成分を直接評価します。

$$
q_{b,j} = \sqrt{s g d^3}\;\Phi\left(\tau_{\ast j}\right)\left[\frac{u_j}{\lvert\mathbf{u}\rvert_j} - \gamma_j\,\frac{z_{b,2(j)} - z_{b,1(j)}}{\delta_j}\right]
$$

$$
\tau_{\ast j} = \frac{g n^2 \lvert\mathbf{u}\rvert_j^2}{H_j^{1/3}\, s g d}, \qquad \gamma_j = \sqrt{\frac{\tau_{\ast c}}{\mu_s \mu_k \tau_{\ast j}}}
$$

決め方の理由は [ADR-0011](adr/0011-edge-bedload.md) にあります。

- 流れの値は、流れを進めた後（時刻 $n+1$）の $u_j$、$\eta$ を使います。
- 辺の水深 $H_j$ は4.1節と同じ規則（風上側の水位 − 両側の高い方の河床高）で、$\lvert\mathbf{u}\rvert_j$ は摩擦と同じ値（法線成分 $u_j$ と、両側のセルの Perot 流速の平均の接線成分）です。したがって掃流力は、運動量式の摩擦が表す底面のせん断力と一致します。
- $H_j <$ `dry_depth` の辺と、$\tau_\ast \le \tau_{\ast c}$ の辺では $q_{b,j} = 0$ です（勾配の項も 0）。
- 母点を結ぶ線分が辺に直交するので、法線方向の河床勾配は2点差分で求まります。勾配補正は全方向（流下方向と横断方向）にかけます。
- 周期境界では河床高を縦断勾配 $I_p$ を除いた値で持つので、勾配の項は $(z_{b,2(j)} - z_{b,1(j)})/\delta_j - I_p n_{x,j}$ にします。
- 壁の辺では $q_{b,j} = 0$ とします。上流端と下流端の辺では、内側のセルの水深と辺の流速から同じ式で求めます。勾配には、内側のセルの河床勾配ベクトル（隣のセルとの差の最小二乗）の法線成分を使います（6章）。
- 辺で評価するので、河床の式は中心差分になります。前進オイラーと中心差分が安定する条件は、河床波の速さ $c$ と勾配補正の拡散係数 $D = M q_b\gamma/(1-\lambda)$ で $c^2\Delta t/(2D) \le 1$ です。`examples/sill.toml` で見積もると、等流部で 0.05、床止めの直下の射流部で 0.5 程度です。

**非侵食層（床止め）**：セル $i$ から1ステップで流出する土砂量 $O_i$ と、非侵食層より上にある土砂量 $S_i$ を比べます。

$$
O_i = M \Delta t \sum_{j\in\partial i,\ s_{i,j} q_{b,j} > 0} s_{i,j}\, L_j\, q_{b,j}, \qquad S_i = (1-\lambda)\, A_i \max\left(z_{b,i} - z_{ne,i},\ 0\right)
$$

$O_i > S_i$ のときは、セル $i$ から流出する辺の流砂量をすべて $S_i/O_i$ 倍に縮めます。各辺の流砂は、上流側のちょうど1つのセルにとってだけ流出になるので、縮める倍率は辺ごとに一意に決まります。

- $S_i$ は時刻 $n$ に積もっている量だけで判定し、同じステップの流入は数えません。反復なしで $z_b \ge z_{ne}$ が保証されます。床止めの上に来た土砂は、1ステップ遅れて通過します。
- セルのループで倍率を求め、辺のループで、流出側のセルの倍率を掛けます。

### 4.7 浮遊砂

浮遊砂の移流には、連続式と同じ θ 平均の単位幅流量 $q_j^{\theta} = H_j^n\left[\theta\, u_j^{n+1} + (1-\theta)\, u_j^n\right]$ を使います。

$$
A_i\,(hC)_i^{\ast} = A_i\,(hC)_i^{n} - \Delta t\sum_{j\in\partial i} s_{i,j}\, L_j\, q_j^{\theta}\, C_{\mathrm{up}(j)}^n + \Delta t\sum_{j\in\partial i} L_j\,\frac{(h\varepsilon)_j}{\delta_j}\left(C_{k(i,j)}^n - C_i^n\right) + \Delta t\, A_i\, E_i
$$

沈降は陰的に扱い、$C_i^{n+1} = (hC)_i^{n+1} / h_i^{n+1}$ を使って次のように解きます。

$$
(hC)_i^{n+1} = \frac{(hC)_i^{\ast}}{1 + \Delta t\, w_f\, r_i / h_i^{n+1}}
$$

- 連続式と同じ流量を使うと、濃度が一様な状態（$C$ = 一定）が移流で崩れません。上の式に $C$ = 一定を代入すると、移流の項が連続式の体積変化と打ち消し合います。
- 水深が `min_depth` 未満のセルでは $E_i = 0$ として浮上を止めます。浅いセルで濃度が過大になるのを防ぐためです。
- 移流は陽的な風上差分なので、$hC \ge 0$ を保つには7章の Courant 条件が必要です。

### 4.8 河床の更新

$$
(1-\lambda)\, A_i\left(z_{b,i}^{n+1} - z_{b,i}^{n}\right) = -M\Delta t\left[\sum_{j\in\partial i} s_{i,j}\, L_j\, q_{b,j} + A_i\left(E_i - w_f\, r_i\, C_i^{n+1}\right)\right]
$$

更新の後、次の処理を行います。

浮遊砂を実装するまでは、浮遊砂の項（$E_i - w_f r_i C_i^{n+1}$）は 0 です。丸め誤差で $z_{ne}$ をわずかに下回ったセルは、$z_{ne}$ にそろえます。決め方の理由は [ADR-0012](adr/0012-bed-update.md) にあります。

1. **安息角補正**：$\lvert z_{b,2(j)} - z_{b,1(j)}\rvert / \delta_j > \tan\phi$ となる辺で、その辺だけならちょうど $\tan\phi$ になる河床の体積の半分（緩和係数 $\omega = 0.5$）を、高い側から低い側へ移します。

$$
V_j = \omega\,\frac{\lvert z_{b,1(j)} - z_{b,2(j)}\rvert - \delta_j\tan\phi}{1/A_{1(j)} + 1/A_{2(j)}}
$$

   - 全辺を同時に評価する Jacobi 型です。辺のループで移す量を辺の配列に書き、セルのループで集めます。両側のセルの値が丸め誤差まで一致するので、体積は厳密に保存されます。
   - 周期境界では、高低差に除いた勾配の分 $I_p n_{x,j}\delta_j$ を戻した、実際の高低差で判定します。
   - 高い側のセルから出る量の合計は、非侵食層より上の量で頭打ちにします。床止めの下流面のように、高い側に土砂のない段差は崩れません。
   - 1回の反復で動いた厚さの最大値が粒径の $10^{-3}$ 倍未満になるまで（最大 100 回）繰り返します。
2. **水位の調整**：水深が `dry_depth` 以上のセルでは水位 $\eta_i$ を保ち、水深が河床の変化分だけ変わります（河床が水面を超えたら $\eta_i = z_{b,i}$）。水深が `dry_depth` 未満のセルでは、水深を保ちます（ほぼ乾いたセルで水を作らないため）。水の体積の変化 $\sum_i A_i (h_i^{\text{新}} - h_i^{\text{旧}})$ を質量収支に記録します（`MassBalance::bed_exchange`）。

## 5. 1ステップの手順

| 段階 | 処理 | 主なモジュール |
|---|---|---|
| 時間刻み 0 | 固定の $\Delta t$、またはセルの流出量の比から決める可変の $\Delta t$（7章） | `model::Model::step`, `hydro::HydroSolver::max_outflow_rate` |
| 流れ 1 | 辺の役割（内部・壁・境界）、辺の水深 $H_j^n$ と wet/dry の判定、辺の単位幅流量 $q_j^n$ | `hydro::edge_depths` |
| 流れ 2 | 移流を陽的に評価して $F u_j^n$ を求める | `hydro::advection` |
| 流れ 3 | 摩擦係数 $\gamma_j$、$G_j$、重み $w_j$、右辺 $b_i$ を組み立てる | `hydro::HydroSolver::step` |
| 流れ 4 | Newton 法＋行列フリー PCG で $\boldsymbol{\eta}^{n+1}$ を求める | `hydro::WaterLevelOperator`, `linsolve::pcg` |
| 流れ 5 | $u_j^{n+1}$ を更新し、θ 平均流量 $q_j^{\theta}$ を保存する | `hydro::HydroSolver::step` |
| 流れ 6 | 乾いたセルの水位を河床高に揃え、Perot 法でセル中心の流速 $\mathbf{u}_i^{n+1}$ を再構成 | `hydro::reconstruct_cell_velocity` |
| 流砂 7 | 辺の掃流砂量 $q_{b,j}$（勾配補正、非侵食層の頭打ち） | `sediment::SedimentSolver::compute_bedload` |
| 流砂 8 | 浮遊砂 $(hC)^{n+1}$（沈降は陰的） | `sediment::SedimentSolver::step_suspended` |
| 流砂 9 | Exner 式（MORFAC）、安息角補正、水位の調整 | `sediment::SedimentSolver::update_bed` |
| 収支 10 | 水と土砂の質量収支を更新し、出力時刻なら VTK を書く | `model`, `io` |

セル中心の流速 $\mathbf{u}_i^n$（摩擦と移流に使う）は、前のステップの最後（初期状態ではモデルの組み立て時）に再構成した値です。流砂の段階（7〜9）は、時刻が `time.morphology_start` を過ぎてから行います。

## 6. 境界条件

| 境界 | 流れ | 流砂 |
|---|---|---|
| 壁 | $u_j = 0$。辺を水位方程式から除く | 掃流砂・拡散フラックスとも 0 |
| 上流端（流量） | 流量 $Q$ を境界辺に通水能 $K_j = L_j H_j^{5/3}/n$ に比例して配分し、既知の流入 $Q_i$ として右辺に入れる。$H_j = \max(\eta_s - z_{b,i},\ $`dry_depth`$)$ で、$\eta_s$ は上流端の湿ったセルの水位の平均（断面で共通の水位） | 平衡給砂（境界辺の流れから求めた掃流砂量を流入させる。下の注）、一定の給砂量（河床変動を始めたときの平衡給砂量の合計を保ち、各ステップの平衡給砂量に比例して辺に配る。上流端の流れが速くなる条件では給砂不足で上流の区間が下がり、流路が分化しうる。[ADR-0013](adr/0013-constant-sediment-supply.md)）、または給砂なし |
| 下流端（水位固定） | 与えた水位 $\eta_b$ を境界辺の上の値とし、母点から辺までの距離 $\delta_j/2$ で水位勾配を評価する。$w_j\eta_b$ を右辺へ移す | 境界辺の流れから求めた掃流砂量で自由に流出（流入はさせない） |
| 下流端（水深一定） | 境界セルの河床高に一定の水深を加えた値を $\eta_b$ とする | 移動床とは組み合わせない（設定の誤りにする。下の注） |
| 周期境界（直線水路の上流端と下流端） | 両端を継いだ辺を内部の辺として扱う。河床高と水位は縦断勾配 $I_p$ を除いた値で持ち、勾配は運動量式の体積力 $g I_p n_x$ にする（下の注） | 継ぎ目の辺も内部の辺と同じ式。給砂・流出はなく、土砂の量は保たれる |

開境界の決め方の理由は [ADR-0007](adr/0007-open-boundaries.md) にあります。

**開境界の掃流砂**（[ADR-0011](adr/0011-edge-bedload.md)）：上流端と下流端の辺の掃流砂量は、内部の辺と同じ式（4.6節）で、水深に内側のセルの水深、流速に辺の流速 $u_j$ を使って求めます。勾配 $\partial z_b/\partial n$ には、内側のセルの河床勾配ベクトル（隣のセルとの差の最小二乗）の法線成分を使います。一様な勾配の等流では、境界辺の流砂量も内部の辺と一致し、平衡河床が保たれます（検証ケース V12）。勾配の項を省くと、境界のセルが系統的に掘れたり堆積したりします。

**下流端の水深一定の境界と移動床**：水深一定の境界では、境界辺の水位が下流端のセルの河床とともに上下します。河床が下がると水面勾配が増えてさらに掘れ、上がると堰上げでさらに堆積するので、正のフィードバックになります。等流に平衡給砂を与えた計算（正方形配置 2 m、地形の時間 30000 s）では、水位固定なら河床は変わりませんが、水深一定では下流端のセルが 600〜900 s の間に 0.93 m 下がり、その後も低下が上流へ進みました。そのため、移動床（`[sediment]`）では水深一定の境界を設定の誤りにし、水位固定を使います（[ADR-0011](adr/0011-edge-bedload.md)）。

- 水位を境界辺の上に置くと、等流が離散解として厳密に成り立ちます（検証ケース V06）。鏡映点に置くと、水面に勾配がある流れで $O(I\delta_j)$ の誤差が出ます。
- 通水能の水深を各セルの水深から求めると、深くなったセルほど流入が増える正のフィードバックで、横断方向に不安定になります。断面で共通の水位を使うと、配分は河床高の違いだけで決まります。
- 流量境界の辺の流速 $u_j$ は、単位幅流量 $Q_j/L_j$ を内側のセルの新しい水深で割った値にします（Perot 再構成と出力に使う）。
- 下流端の辺の水深の風上は、流出なら内側のセルの水位、流入なら $\eta_b$ で決め、河床高には内側のセルの値を使います。

**周期境界**（[ADR-0014](adr/0014-periodic-channel.md)）：交互砂州の検証（V09）のために、直線水路の上流端と下流端をつなぎます（`[boundary.upstream]` と `[boundary.downstream]` をともに `kind = "periodic"` にする）。

- 格子は `mesh::channel::build_periodic` で作ります。六角形配置は x 方向に周期的な三角格子にして、継ぎ目のセルも内部と同じ六角形にします（継ぎ目はジグザグ）。継ぎ目の辺には、反対側の母点の位置への平行移動（`edge_shift`）を持たせ、Perot 法・河床勾配の最小二乗・直交性の検査では、その位置を使います。
- 河床高と水位は、縦断勾配 $I_p$（`[bed] slope`）を除いた値で持ちます。実際の値は、それぞれに $I_p(L - x)$ を足したものです。辺の水深の風上と高い方の河床は、除いた勾配を戻した実際の高さで選ぶので（4.1節、[ADR-0015](adr/0015-periodic-edge-depth.md)）、同じ状態の辺の水深は勾配のある開境界の水路と一致します。
- 給砂と流出はなく、水と土砂の量は保たれます。流量は水の量と摩擦で決まり、河床が変わると変わります。`discharge` は、`Model::new` が初期の等流（水深 $h_n = (nq/\sqrt{I_p})^{3/5}$）を作るための流量です。
- 周期 L の水路と、同じ初期状態を2回並べた周期 2L の水路の解は、丸め誤差の範囲で一致します（`tests/periodic_channel.rs`）。

## 7. 時間刻みの制約

- **重力波**：水位勾配を陰的に扱うので、重力波の Courant 条件 $(\lvert\mathbf{u}\rvert + \sqrt{gh})\Delta t / \Delta x \le 1$ には制約されません。
- **移流**：移流は陽的な風上差分なので、セルから1ステップで流出する体積がセルの体積を超えないことが必要です。おおよそ $\lvert\mathbf{u}\rvert\Delta t/\Delta x \lesssim 1$ です。
  - 運動量の移流がなくても、この制約はかかります。連続式の辺の水深 $H_j^n$ を時刻 $n$ の風上側からとるので、水深を陽的な風上差分で移流しているのと同じになるためです。
  - この条件が成り立てば、運動量の移流の局所クーラン数 $\sigma_j$（4.4節）は 1 以下です。`StepReport` の `max_courant` と `limited_edges` で確かめられます。
  - 段波（速さ $\lvert u\rvert + \sqrt{gh}$）が1ステップで1セル以上進む時間刻みでも安定ですが、時間について1次の誤差が出ます。V06 の初期の過渡状態（段波のクーラン数 1.2）では、段波の先端の流速が約 2 割大きく出ました（[ADR-0008](adr/0008-advection-discretization.md)）。
  - V06（流速 1.87 m/s、格子幅 2 m）では、クーラン数 0.47 で安定、1.9 で不安定でした。0.93 では、正方形配置は安定、六角形配置は不安定でした（六角形のセルは体積に対する流出量が大きく、上の条件を満たさないセルがある）。

$$
\Delta t \sum_{j\in\partial i,\ s_{i,j} q_j > 0} s_{i,j}\, L_j\, q_j \le V_i
$$

- **河床**：MORFAC を掛けた河床の変化が、1ステップで隣のセルとの高低差を逆転させない程度に抑えます。時間刻みは流れのクーラン数だけで決め、1ステップの河床の変化の最大値を記録して表示します（`BedUpdateReport::max_change`）。

### 時間刻みの決め方

設定の `[time]` に、`dt`（固定）か `courant` と `dt_max`（可変）のどちらか一方を書きます。決め方の理由は [ADR-0010](adr/0010-adaptive-time-step.md) にあります。

- **固定**：$\Delta t$ = `dt`。時刻は「開始時刻 + ステップ数 × $\Delta t$」で計算し、足し算の誤差がたまらないようにします。
- **可変**：各ステップの初めに、時刻 $n$ の流量で、上の条件の左辺と右辺の比（セルの流出量の比）$r$ を求め、$\Delta t = \min(C/r,\ \Delta t_{\max})$ とします（流出がなければ $\Delta t_{\max}$）。

$$
r = \max_i \frac{1}{V_i}\sum_{j\in\partial i,\ s_{i,j} q_j^n > 0} s_{i,j}\, L_j\, q_j^n
$$

  - $C$ = `courant` は目標のクーラン数（$0 < C \le 1$）、$\Delta t_{\max}$ = `dt_max` は上限です。上限は、流れのない状態（乾いた状態から始めるとき）で $\Delta t$ を決めるために必須です。
  - 次の区切りの時刻（出力時刻、`morphology_start`、終了時刻）までの残りが $\Delta t$ 以下なら、ちょうどそこまで進めます。残りが $2\Delta t$ 未満なら、残りの半分だけ進めます。
  - 六角形配置のセルは、体積に対する流出量が大きいので、同じ流速でも $r$ が大きくなります。`examples/sill.toml` の定常流（床止めの直下の射流 約 5.7 m/s、格子幅 1 m）では、正方形配置で $r$ = 5.7 1/s、六角形配置で 8.1 1/s でした。$C$ = 0.5 のとき、六角形配置で $\Delta t$ ≈ 0.062 s です。
- 実際のクーラン数 $\Delta t\, r$ は、固定のときも含めて `StepReport::outflow_courant` に記録します。

半陰解法で得をする度合いは、フルード数が小さいほど大きくなります。たとえば `examples/sill.toml` の等流部（$Fr \approx 0.58$）では、重力波と移流の Courant 条件の差は約 2.7 倍にとどまります。床止めの直下の射流（$Fr \approx 3$）では、移流の条件のほうが厳しくなります。この例題では、時間刻みの利点よりも wet/dry と急変流を安定に扱えることが主な目的です。

## 8. 未決事項

実装の過程で決め、結果を ADR に記録します。

- 移流の2次精度化と制限関数（重みは [ADR-0008](adr/0008-advection-discretization.md)、形は [ADR-0009](adr/0009-advection-face-velocity.md) で決定）
- 運動量保存形とエネルギー水頭保存形を切り替える判定方法
- セル内の河床の起伏（サブグリッド地形、Casulli, 2009）を体積関数に入れるかどうか
- 二次流による掃流砂の向きの補正（湾曲水路の検証ケース V10 で必要になる）
- 混合粒径（平野, 1971 の交換層モデル）
- 水平渦粘性の要否

## 参考文献

書誌情報と入手方法は [docs/references.md](../references.md) にまとめています。
