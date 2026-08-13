# Reconciled speciation contrast の数式と最小例

この文書は、NWKIT が reconciled speciation contrast をどのように作るかを、
実装に対応した数式と手計算できる例で示す。コマンド全体の使い方、入力形式、
進化モデル、replicate、PGLS の選択肢については
[`RECONCILED_CONTRASTS.md`](RECONCILED_CONTRASTS.md) を参照すること。

## 一文での定義

Reconciled speciation contrast は、**遺伝子木上の通常の phylogenetic
independent contrast (PIC) のうち、reconciliation によって、ある種木の
種分化イベントの異なる娘クレードを比較すると確認されたもの**である。

Reconciliation 自体が発現量を差し引くわけではない。役割分担は次のとおりで
ある。

1. reconciliation が、遺伝子木ノードを種木ノードへ写像する。
2. reconciliation が、採用可能な種分化ノードと numerator/denominator の向きを
   決める。
3. PIC recursion が、遺伝子木の枝長を使って発現量の contrast を計算する。
4. 同じ種木ノードで、種木の枝長を使って説明変数の contrast を計算する。
5. `species_event_id` と娘クレードIDで両者を結合し、PGLSを当てる。

したがって、応答contrastと説明変数contrastは同じ種分化を表すが、祖先値と
contrast分散はそれぞれ異なる木上で計算される。

## 1. Reconciliation によるイベント選択と符号の固定

遺伝子木の内部ノードを $g$、そこに写像された種木ノードを $s$ とする。
$s$ の2つの娘クレードが含む種集合を $S_N,S_D$ と書く。NWKIT は、各クレードの
ソート済み子孫種名を辞書順に並べ、先の側を numerator ($N$)、後の側を
denominator ($D$) とする。この規則は Newick 内の子ノード順に依存しない。

$g$ の2つの直下の子が含む種集合を $G_1,G_2$ とする。次のどちらかが成り立つ
ときだけ、$g$ は contrast に使える種分化イベントとなる。

$$
(G_1\subseteq S_N \land G_2\subseteq S_D)
\quad\text{or}\quad
(G_2\subseteq S_N \land G_1\subseteq S_D).
$$

該当する遺伝子木の子をそれぞれ $g_N,g_D$、それらが含む種集合を
$G_N,G_D$ と向き付ける。ノード自体がresolved speciationでない場合、両方の
遺伝子木の子が同じ種木娘クレードに入る場合、または写像が曖昧な場合は
`eligible=no` となり、reconciled speciation PGLS には使われない。

完全coverageはさらに

$$
G_N=S_N,\qquad G_D=S_D
$$

を要求する。部分集合ではあるが等しくない場合は、イベントとしては解釈可能でも
`coverage_status=partial` となる。既定の `--speciation-coverage complete` はこれを
除外し、`any` は明示的な感度分析として含める。

種木ノード $s$ の子孫種集合から作られる安定IDが `species_event_id` である。
$S_N,S_D$ のIDも `species_numerator_event_id` と
`species_denominator_event_id` として保存される。`species_event_id` が遺伝子木と
種木のcontrastの結合キーになり、2つの娘クレードIDが向きの一致を検証する。

## 2. PIC recursion

以下の式は発現量にも種形質にも共通である。形質を $z$ とし、選択した進化モデル
で変換した枝 $e$ の分散長を $b_e$ とする。Brownian model なら $b_e$ は入力枝長
そのものである。kappa、delta、OUなどでは、まず等価な $b_e$ に変換してから
同じrecursionを使う。

### Tip

tip $i$ では、推定値と下から持ち上げる分散を

$$
\hat z_i=z_i,\qquad u_i=b_i
$$

と置く。$b_i$ はtipへ入る枝の変換後分散長である。

### Internal node

内部ノード $g$ の向き付けられた子を $g_N,g_D$ とし、それまでに計算された値を
$(\hat z_N,u_N)$ と $(\hat z_D,u_D)$ とする。NWKIT が出力する4つの主要量は

$$
\begin{aligned}
c_g &= \hat z_N-\hat z_D
&&\text{(`raw_contrast`)},\\
q_g &= u_N+u_D
&&\text{(`contrast_variance`)},\\
c_g^{\mathrm{std}} &= \frac{c_g}{\sqrt{q_g}}
&&\text{(`standardized_contrast`)},\\
\hat z_g &=
\frac{u_D\hat z_N+u_N\hat z_D}{u_N+u_D}
&&\text{(`ancestral_estimate`)}.
\end{aligned}
$$

$\hat z_g$ は分散の逆数で重み付けした祖先値である。さらに上位ノードへ渡す分散は

$$
u_g=b_g+\frac{u_Nu_D}{u_N+u_D},
$$

となる。$b_g$ は $g$ から親へ入る枝の変換後分散長であり、rootでは0である。
この祖先値と調整済み分散をpostorderで繰り返すため、深い内部ノードのcontrastは
単純なtip平均の差とは限らない。

### 線形係数としての表現

tip値をベクトル $\mathbf z$ とすると、各raw contrastは

$$
c_g=\boldsymbol\ell_g^\mathsf T\mathbf z
$$

と書ける。$\boldsymbol\ell_g$ はrecursionから得られるPIC係数である。選択した
$r$ 個の遺伝子木contrastを行に並べた行列を $L$ とすると

$$
\mathbf c=L\mathbf z.
$$

tip平均にsampling covariance $D$ がある場合、biological replicate由来の
contrast covarianceは

$$
M=LDL^\mathsf T
$$

である。NWKITは対角成分だけでなく、共有tipやbatch補正によって生じる非対角成分
も保持する。tip誤差が独立で
$D=\operatorname{diag}(d_1,\ldots,d_n)$ のときは

$$
U=L\operatorname{diag}(\sqrt{d_1},\ldots,\sqrt{d_n}),
\qquad M=UU^\mathsf T
$$

と因子のまま保持できる。NWKITはcontrast数が500を超える場合、この $U$ の
非ゼロ要素を `covariance_representation=factor-loading` として出力し、PGLSでも
$M$ を密行列化せずに用いる。batch補正がある場合も、tip covarianceを

$$
D=\operatorname{diag}(\mathbf d)+HH^\mathsf T
$$

という「対角＋低ランク」形式で直接計算し、contrast側では
$U=[L\operatorname{diag}(\sqrt{\mathbf d}),LH]$ とする。出力時には各行の最大loading
の16 machine epsilon以下だけを除く。この丸め誤差水準の疎化を除けば上式と同じ
共分散であり、tip数の二乗の行列は作らない。

## 3. 遺伝子contrastと種形質contrastの対応

遺伝子木の採用ノードを $g_i$、対応する種木イベントを $e(i)$ とする。
応答である発現量のraw contrastを

$$
y_i=c^{\mathrm{gene}}_{g_i}
$$

とする。説明変数 $j$ の種木raw contrastは、種木イベント $e$ で一度だけ計算して

$$
x_{e,j}=c^{\mathrm{species}}_{e,j}
$$

とする。PGLSの行 $i$ には $x_{e(i),j}$ が割り当てられる。

重要なのは、最終的なreconciled PGLSが `standardized_contrast` 同士を回帰する
のではない点である。raw contrastを使い、遺伝子PICの進化分散を

$$
G=\operatorname{diag}(q_1,\ldots,q_r)
$$

としてcovariance modelへ入れる。したがって、基本の固定効果モデルは切片なしで

$$
\mathbf y=X\boldsymbol\beta+\boldsymbol\varepsilon
$$

であり、$X_{ij}=x_{e(i),j}$ である。既定のhierarchical modelでは概略

$$
\operatorname{Cov}(\boldsymbol\varepsilon)
=\sigma^2G+M
+\tau_{\mathrm{event}}^2Z_{\mathrm{event}}Z_{\mathrm{event}}^\mathsf T
+\text{lineage random-slope component}
$$

を推定する。

### 3.1 説明変数にもreplicateがある場合

説明変数 $j$ の観測された種木contrastを $\hat{\mathbf x}_j$、真の潜在contrastを
$\mathbf x_j^*$、replicateからPIC伝播したsampling covarianceを $M_{x,j}$ とする。
種木contrastの進化分散を対角に並べた $K_{x,j}$ と進化rate $\rho_j$ を使い

$$
\hat{\mathbf x}_j=\mathbf x_j^*+\boldsymbol\eta_j,
\qquad
\boldsymbol\eta_j\sim N(\mathbf0,M_{x,j}),
$$

$$
\mathbf x_j^*\sim N(\mathbf0,\rho_jK_{x,j})
$$

と置く。$\rho_j$ は説明変数だけのmarginal likelihoodで推定する。Gaussian
conditioningにより

$$
\mathbf m_j=E(\mathbf x_j^*\mid\hat{\mathbf x}_j)
=K_j(K_j+M_{x,j})^{-1}\hat{\mathbf x}_j,
$$

$$
S_j=\operatorname{Var}(\mathbf x_j^*\mid\hat{\mathbf x}_j)
=K_j-K_j(K_j+M_{x,j})^{-1}K_j,
\qquad K_j=\rho_jK_{x,j}.
$$

$R$ を「一つの種木イベント」を「それに対応する複数の遺伝子木行」へ写す行列と
すると、固定効果の説明変数は $R\mathbf m_j$、条件付き応答分布は

$$
\mathbf y\mid\hat X
\sim N\left(\sum_j\beta_jR\mathbf m_j,
V+\sum_j\beta_j^2RS_jR^\mathsf T\right)
$$

となる。$V$ は上記の応答側covarianceである。$\beta_j^2RS_jR^\mathsf T$ により、
説明変数誤差を応答誤差として単純加算せず、回帰係数に依存する形で正しく伝播する。
同じ種分化を通るパラログは同一の潜在値と共分散を共有するため、ここでも
pseudoreplicateにはならない。説明変数にsampling covarianceがない場合は
$M_{x,j}=0,\ \mathbf m_j=\hat{\mathbf x}_j,\ S_j=0$ となり、従来の式へ戻る。

この条件付きcovarianceは $\beta_j$ に依存するため、標準的なREMLの導出条件を
満たさない。したがって説明変数誤差を含むモデルはMLで当てはめる。カテゴリカル
説明変数のbiological replicatesについては、factor coding後の一観測の共分散を
$\Omega_i$、独立なreplicate数を $n_i$ として、平均coding値の共分散
$\Omega_i/n_i$ を全cross-column項ごと伝播する。これは離散状態の全組合せを
厳密に周辺化するものではなく、linear predictor上のGaussian moment近似である。

## 4. 同じ種分化を通る複数パラログ

古い重複の後に同じ種分化を複数のパラログ系統が通ると、異なる $g_i$ が同じ
$e(i)$ を持つ。このとき $X$ の同じ種形質contrastが複数行に現れるが、独立な
種分化の反復が増えたわけではない。

イベント $e$ に属する遺伝子contrast数を $k_e$ とすると、既定の
`--event-weighting equal` は、対角だけの単純な場合には各行へ $1/k_e$ の情報重み
を与えることに相当する。実装上は $B_{ii}=\sqrt{k_{e(i)}}$ としてworking
covarianceを

$$
V_{\mathrm{work}}=BVB
$$

へ変換する。Gaussian ML/REMLでは行数そのものを尤度の標本数にせず、イベント数
$m$ を有効標本数とする。$n$ を遺伝子contrast行数として、共分散のlog determinant
項は

$$
\frac{m}{n}\left\{\log|V_{\mathrm{work}}|-2\log|B|\right\}
$$

とする。$-2\log|B|$ は既知の行スケーリングのJacobianを相殺する項である。したがって
同一のparalog行を何コピー追加しても、`replicate-reml` の係数、進化rate、標準誤差、
log likelihoodは変わらない。parametric bootstrapとlineage likelihood-ratio検定も
同じevent-level擬似尤度で再fitする。

さらに、識別可能なら同じ `species_event_id` を共有する行に
species-event random effectを入れる。残差自由度も遺伝子contrast数ではなく

$$
\mathrm{df}=n_{\mathrm{species\ events}}-p
$$

で計算する。ここで $p$ は固定効果の数である。

## 5. 最小の手計算例：2種、2パラログ

種木、遺伝子木、値を次のようにする。枝長はすべて1、進化モデルはBrownianとする。

```text
species tree:  (A:1,B:1);

gene tree:     ((A_g1:1,B_g1:1):1,
                (A_g2:1,B_g2:1):1);
                 \____ copy 1 ____/  \____ copy 2 ____/
                    speciation           speciation
                         \______________/
                            duplication
```

種形質と発現量は次のとおりとする。

| tip | 値 |
|---|---:|
| species trait $x_A$ | 2 |
| species trait $x_B$ | 6 |
| expression $y_{A,g1}$ | 3 |
| expression $y_{B,g1}$ | 7 |
| expression $y_{A,g2}$ | 10 |
| expression $y_{B,g2}$ | 14 |

種木の娘クレードは辞書順で $A$ がnumerator、$B$ がdenominatorになる。
2つの遺伝子木speciation nodeも同じ向きに固定され、同じ
`species_event_id`（ここでは $E_{AB}$ と略記）へ写像される。rootはduplication
なのでspeciation contrastから除外される。

### 種形質側

$$
\begin{aligned}
c_x &= x_A-x_B=2-6=-4,\\
q_x &= 1+1=2,\\
c_x^{\mathrm{std}} &= -4/\sqrt{2}=-2\sqrt{2}\approx-2.828,\\
\hat x_{AB} &= (2+6)/2=4.
\end{aligned}
$$

### 遺伝子発現側

copy 1では

$$
c_{y,1}=3-7=-4,\qquad q_{y,1}=2,
\qquad c_{y,1}^{\mathrm{std}}\approx-2.828.
$$

copy 2でも

$$
c_{y,2}=10-14=-4,\qquad q_{y,2}=2,
\qquad c_{y,2}^{\mathrm{std}}\approx-2.828.
$$

PGLSへ渡る最小表を、人が読めるIDに略記すると次のようになる。

| gene node | `species_event_id` | response `raw_contrast` | response `contrast_variance` | predictor `raw_contrast` |
|---|---|---:|---:|---:|
| copy 1 speciation | $E_{AB}$ | -4 | 2 | -4 |
| copy 2 speciation | $E_{AB}$ | -4 | 2 | -4 |

形式的には両行とも $-4=\beta(-4)$ なので $\hat\beta=1$ である。しかし、これは
独立な種分化が1個しかない例である。2行を独立な $n=2$ と数えてはいけない。
このデータ集合のカウントは、fit成立時の次の診断列に対応する。

```text
n_gene_contrasts = 2
n_species_events = 1
n_repeated_gene_contrasts = 1
```

となる。また説明変数が1個なら $\mathrm{df}=1-1=0$ なので、NWKITはこの例だけで
PGLS推論を行うことを拒否する。この例はcontrast演算を示すための最小例であり、
実際の推論には、説明変数の数より多い**独立な種分化イベント**が必要である。

## 6. 最小例を低レベルコマンドで再現する

上の例は次の5ファイルで再現できる。

`species.nwk`

```text
(A:1,B:1);
```

`gene.nwk`

```text
((A_g1:1,B_g1:1):1,(A_g2:1,B_g2:1):1);
```

`gene-to-species.tsv`

```tsv
leaf_name	species_label
A_g1	A
B_g1	B
A_g2	A
B_g2	B
```

`expression.tsv`

```tsv
leaf_name	expression
A_g1	3
B_g1	7
A_g2	10
B_g2	14
```

`species-trait.tsv`

```tsv
leaf_name	body_size
A	2
B	6
```

まずreconciliationを作り、次に両方のcontrastを計算する。

```sh
nwkit reconcile \
  --infile gene.nwk \
  --species-tree species.nwk \
  --species-map-tsv gene-to-species.tsv \
  --tree-id OGmini \
  --outfile reconciliation.tsv

nwkit contrast \
  --infile gene.nwk \
  --trait expression.tsv \
  --columns expression \
  --reconciliation reconciliation.tsv \
  --event-type speciation \
  --tree-id OGmini \
  --outfile gene-contrasts.tsv

nwkit contrast \
  --infile species.nwk \
  --trait species-trait.tsv \
  --columns body_size \
  --tree-id species \
  --outfile species-contrasts.tsv
```

`gene-contrasts.tsv` には2行のeligible speciation contrastが残り、両方の
`species_event_id` が一致する。`species-contrasts.tsv` の対応行では
`branch_clade_id` がそのIDと一致し、numerator/denominator clade IDも一致する。
実際のIDは子孫tip名から作る `clade-sha256:...` であり、Newickの子順を変えても
変わらない。

## 7. 出力を読むときのチェックリスト

- `event_type=speciation`, `mapping_status=mapped`, `eligible=yes` か。
- 主解析なら `coverage_status=complete` か。
- gene側の `species_event_id` とspecies側の `branch_clade_id` が一致するか。
- numerator/denominator のspecies clade IDも一致するか。
- 回帰に使う値は `raw_contrast`、進化分散はgene側の
  `contrast_variance` であることを取り違えていないか。
- `n_gene_contrasts` ではなく `n_species_events` を独立な比較数として見ているか。
- 同一イベントのパラログが多い場合、`event-weighting=equal` とevent random
  effectの診断を確認したか。

この構成により、パラログごとの発現差を潰さずに保持しつつ、同じ種分化contrastを
繰り返し割り当てたことによるpseudoreplicationを最終モデルで明示的に扱える。

## 8. 種分化後に形質を獲得し、片方のパラログだけが応答する場合

固定効果だけの式 $y_i=x_{e(i)}\beta+\epsilon_i$ は、同じ種分化イベントを通る
すべてのパラログに共通の平均傾き $\beta$ を推定する。NWKITのhierarchical modelは
これを

$$
y_i=x_{e(i)}\{\beta+b_{\ell(i)}\}+a_{e(i)}+\epsilon_i,
\qquad b_\ell\sim N(0,\tau_{\mathrm{lineage}}^2)
$$

へ拡張する。$b_\ell$ は `lineage_clade_id` ごとの傾き偏差、$a_e$ は同じ
`species_event_id` に共有される応答偏差である。したがって、古い重複で生じた
copy 1だけが形質と同方向に発現変動しても、共通の $\beta$ へ完全に平均化せず、
copy 1のtotal slope $\beta+b_1$ とcopy 2のtotal slope $\beta+b_2$ を別々に出せる。

たとえば4個の独立な種イベントで $x=(-1,-1,1,1)$、copy 1の発現contrastが
$y_1=(-2,-2,2,2)$、copy 2が $y_2=(0,0,0,0)$ なら、概念的には

$$
\hat\beta\approx1,\qquad
\hat b_1\approx+1,\qquad
\hat b_2\approx-1,
$$

ゆえにtotal slopeはcopy 1で約2、copy 2で約0となる。実際の推定値は進化分散、
sampling covariance、event weightingと部分プールにより0へ縮小される。
`random-effects.tsv` の `conditional_mode` が $\hat b_\ell$、
`total_coefficient` が $\hat\beta+\hat b_\ell$、`reliability` がデータから
得た系統別情報の強さを表す。区間は推定済み分散成分を条件とする
empirical-Bayes区間なので、分散成分そのものの不確実性までは含めない。

平均効果が0でも正負の系統別応答が相殺している可能性がある。このため
`--lineage-inference` は次の2帰無仮説を分ける。

$$
H_{0,\mathrm{het}}:\tau_{\mathrm{lineage}}^2=0
$$

と、説明変数群 $g$ について

$$
H_{0,\mathrm{joint}}:\beta_g=0
\quad\text{かつ}\quad
\text{その群のlineage slopeを持たない}.
$$

前者だけなら分散境界の50:50混合カイ二乗近似を使用できる。後者は通常の
カイ二乗分布にならないため、P-valueにはparametric bootstrapを使う。

なお、この拡張でも形質獲得と発現変化の枝内の時間順序は直接推定しない。
カテゴリカル形質についてはraw-input modeのMk stochastic mappingが各枝の
`from_state -> to_state` posterior frequencyを推定し、その枝より下の種分化
イベントをまとめて除くorigin leave-one-outを提供する。これは
「その起源に依存して結論が出ているか」を調べる診断であり、発現変化の時刻を
同定する因果モデルではない。
