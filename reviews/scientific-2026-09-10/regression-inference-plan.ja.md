# 回帰推論の対応計画（レビュー第1・4・6項）

2026-09-10。計画作成後の「実装して」に基づき、割当worktreeで本番コード・テスト・CLI・文書を更新した。実装範囲と検証結果は [実装記録](regression-inference-implementation.ja.md) を参照。原本レビューと他タスクは変更せず、commit/push・大規模校正実験は実施していない。

## 1. 結論と参照状態

目的関数だけの局所修正では足りない。**生物学的生成モデル、イベント均衡化による推定対象、観測機構、正則化、検定の帰無分布を明示的に分け、それらを一つの推論契約で結ぶ**。イベント平均βを既定の主結果とし、その推定方程式に合う不確実性計算を実装する。通常の階層尤度による共通β・分散・lineage推論は明示的な補助解析とする。打ち切りbootstrapと正則化検定は、この同じ生成・再推定基盤に載せる。

### ユーザー確認済みの方針

推定対象について説明・推奨後、ユーザーから「その方針で進めて」と指示を受けた。以下は確定事項とし、再度の選択確認は不要。

- 既定の主推定対象はイベント平均β_E。種分化イベントごとの総重みを等しくし、パラログ数の多いイベントへの偏りを抑える。
- 共通βの正規化階層モデルは補助解析として明示選択できる。分散・lineageの検定がどのモデルに基づくかを分離して報告する。
- イベント均衡化の意図は維持するが、旧疑似尤度・逆Gram共分散・working covarianceからの生成は維持しない。
- 共通βとの結果差は診断として報告する。結果を見て有意になる方へ自動切替したり、二つの推定量を混ぜて主結果を作ったりしない。
- その後の明示的な「実装して」により本番コード変更へ進めた。commit/pushと大規模実験は実施していない。

- 作業場所：`/Users/kf/.codex/worktrees/39bb/nwkit`。開始時clean、detached HEAD `8d3a4a487bc57b1363c957e1e2bc1e2bb7fab5d6`。AGENTS.md確認済み。branch変更・fetchは不要な読取計画作業として行っていない。
- 原本：`/Users/kf/repos/nwkit`。確認時HEADは同じだが多数の未コミット変更・新規ファイルがある。割当worktreeにレビューがないため、原本の `reviews/scientific-2026-09-10/{REVIEW.md,reproduce.py,results.json}` を読取専用で参照した。
- `regress.py` SHA-256は両場所・レビューとも `b6d77775ba510e90e173b01f2b62bef1203a3d58a28abf95fdcb934c4a1db2cd`、`gaussian.py` は `a2a4eb7e58b811f06c7f9e98a6797f40480c9b23ab7a23b5bac3402faf355eea`。
- `phylogenetic_glmm.py` はworktreeが `3a154cac56be4de5bfeb927b3d79edb7c1db2c034760633b381ffb46c6753d15`、原本がレビューと同じ `45b12f49afd0e2cf3606c32e9e79ac803cbfff7f2bb323c6539eeda312ec92e6`。差分はscalar random modeのNewton polishingとscore判定。担当反例の生成器・検定関数は変わらず、両場所で再現した。実装時はこの既存最適化変更を保全する。
- `measurement_error.py` は両場所 `604a60912b4ce93570442c955c2eadff8145f5370d8979943896987688e2b7b6`、`regression_pipeline.py` は両場所 `1ebc1788f5c696e5b6cdc8df5d565a473f157545ca3ffc38ee60d8db6de23277`。
- 原本の未追跡 `regression_selection.py` 等と変更中のCLI・回帰文書は統合時に再確認する。ここでの行番号はworktreeの参照状態に対応する。

## 2. 確認できた問題と適用範囲

### 2.1 イベント疑似尤度

`regress.py:2271` は行固有進化分散と応答sampling covarianceをイベント内行数で膨らませ、共有event/lineage効果はそのままにする。`:1687` は完全なlog determinantをイベント単位の量に交換し、二次形式はworking covarianceの逆行列のまま。event variance成分の有無によってdeterminantの定義自体が切り替わる。`:1717` はGLSの逆Gram行列を係数共分散とする。同じ疑似determinantは `measurement_error.py:1471` にもある。

独立な20イベント、各2行、平均既知0、真のブロック共分散 `I+11'`、working covariance `2sI+t11'` では、期待負対数目的関数/イベントは定数を除き

`f(s,t)=1/2 { log((2s+t)/2) + 1/(2s) + 3/(2s+2t) }`。

真値(1,1)の20イベント合計勾配 `(-2.083333,-0.416667)`、数値最適値 `(1.3660261,1.0000005)` を再現した。解析的にも `t=1, s=(1+sqrt(3))/2=1.3660254` が停留点となる。これは真の生物学的残差分散を推定する目的関数としての反例であり、全モデルの係数biasや第一種過誤率を測った結果ではない。working scaleの別の推定対象と定義し直せても、現行の生物学的分散としての出力・生成への流用は正当化されない。

別の反例はevent効果なし、10イベント各2行、切片1個、ML。現行rate推定は `RSS/(2m)`、bootstrap生成分散は `2*rate` なので再fit期待比は `(n-p)/m=1.9`。元rate `0.70039256`、期待再fit `1.33074587` を確認。通常の正規化Gaussian MLであれば、同じn=20,p=1の再fit期待比は **0.95** である。修正合格条件を「ML bootstrap平均が厳密に1倍」にしてはいけない。REMLの不偏性は別に確認する。

`regress.py:2050` の半χ²、係数のt基準、bootstrapの中心化tailには現行疑似目的関数の下で校正保証がない。`regression_pipeline.py:1707` のshape再fit用tip再構成にも `fitted_covariance_factor` の流用が及ぶ。通常のcontrast-space bootstrapだけ直すのでは不十分。

なお、正しく指定された平均モデルでは別のworking weightでもβが一致推定になる場合がある。「全ての係数がbiasを持つ」とは結論しない。少数event、交差するlineage、共有測定誤差下の有効独立数と被覆率は未検証。

### 2.2 打ち切りbootstrap

`phylogenetic_glmm.py:2070` は正規乱数を引いた後、元censored行だけをNaNに戻す。生成値でラベルを更新せず、元exact行にも検出制約を適用しない。現在のboundsは観測された区間を記述し、全行の検出装置・観測規則を記述するものではない。

既知SD1、潜在効果分散0の極限、検出限界0、元50 censored/50 exactの例でexact行50個のうち25個が負となった。期待scoreは `-50*phi(0)/Phi(0)=-39.894228`。正しい左打ち切り生成なら1行あたり `Phi(0)*[-phi(0)/Phi(0)] + integral_0^infinity y phi(y)dy=0` で、独立数値積分は `-5.6e-16` だった。

これは現行bootstrapの反例。`_censored_gaussian_log_likelihood` のdensity/CDF/区間確率そのものを否定する証拠ではない。informative censoring、未知の検出限界、採取されない個体のtruncationは別モデルである。

### 2.3 正則化推論

`phylogenetic_glmm.py:1881` は罰則を含むobjective差に通常χ²のtailとprofile閾値を適用する。Waldは主に同objectiveの数値曲率を使う。reported log likelihoodで罰則を引き戻す処理はあるが、罰則付き最適値で評価した無罰則likelihoodは無罰則MLEのlikelihoodと同じではない。

`Y~N(beta,1)`、Gaussian penalty precision3では `beta_hat=Y/4`、差分統計量 `Y²/4~χ²_1/4`。Y=2で現行p=0.3173105、正しいtail=0.0455003を両場所で確認した。罰則曲率の逆数は1/4だが推定量の頻度論的分散は1/16であり、さらに非零βでは縮小biasがある。sandwich分散だけでも被覆は保証されない。

これはGaussian penaltyを渡した汎用関数の反例で、既定Student-t scale2.5の全GLMMが同じ倍率・方向で誤るという主張ではない。Laplace積分誤差、separation、分散境界は独立に残る。

追加確認：現行 `parametric-bootstrap` の係数pはfull fitから生成した係数を中心化してtailを数える（GLMM `:2181`、reconciled `regress.py:3090`）。帰無制約下のpenalized統計量bootstrapと同一ではないため、単にこのoptionへ誘導して完了としてはいけない。

## 3. 採用する統計的定義

### 3.1 共通の生成モデル

樹形・reconciliation・design・既知の測定誤差を条件として、Gaussian応答を

`y = X beta + Z_event a + sum_s Z_s b_s + epsilon + measurement_error`

`C(theta)=sigma² G + M_y + tau_event² Z_event Z_event' + sum_s tau_s² Z_s Z_s'`

とする。イベント内行数kは観測数であり、生物学的分散をk倍する理由にはしない。lineageのwhiteningは単位・factor coding不変性を保ち、元係数単位への逆変換を維持する。

推定対象を二つ区別する。

1. **共通の生物学的関連βと分散成分**：正規化階層likelihoodで推定する。log determinant・二次形式・正規化定数は同じCと行数nを使う。βに依存しないCと固定designの場合に限り正規REMLを選べる。fixed effectの異なるモデル比較はMLで行う。
2. **イベントを等しく重視した平均関連β_E**：母集団として「イベントを一様に選び、そのイベント内の対象paralogを一様に選ぶ」を定義する。モデル不適合・効果異質性がある場合はlikelihoodのβと異なる。完全に同じ観測の複製と、新しい独立生物学的観測の追加は区別する。

係数・分散・random effect・検定に `estimand` と推論basisを記録する。「イベントに同じ情報量」は一般には保証できない。xの大きさ、残差分散、共有効果、測定誤差によって情報量は変わるため、保証するのは事前に定義した損失/score上の総重みとする。

### 3.2 イベント均衡化の選択肢

| 案 | 長所 | 制約と判断 |
|---|---|---|
| A. 正規化階層ML/REML | 分散・係数・生成器・lineage比較を一つの分布で定義できる | 同数イベントへの厳密な総重み制限はない。独立paralogの追加で精度が上がりうる。共通β・分散推論の基盤として採用 |
| B. 明示的イベント均衡推定方程式 | β_Eというイベント平均の目的と行コピー不変性を維持できる | ML/REMLと呼ばない。nuisance推定を含むsandwichまたはモデルbootstrapが必要。共有lineage等のためevent独立cluster-HC1を無条件には使えない。均衡推定の推奨案 |
| C. 正規化されたblock/pair密度のweighted composite likelihood | 各scoreの期待0を導出でき、依存パラメータもpairで識別できる | weightsの設計、交差lineageのpair、Godambe推定が複雑。βの厳密なイベント平均目標と同じとは限らない。Bが不十分な場合の研究代替 |
| D. event平均へ集約し完全な `A C A'` でGLS | 次元が減り、正規化likelihoodを容易に定義できる | event内変動とlineageの情報を失い、分散成分が識別不能になる場合がある。平均関連だけの限定機能 |

Bの最小定義は、事前の単位変換後のraw contrastについて `U_beta=sum_e (1/k_e) sum_i x_ei (y_ei-x_ei' beta)=0`。このlossはコピーした行の重みも再正規化する。分散による追加重みを入れるならβ_Eの定義自体に明記する。nuisance分散はCを用いた不偏Gaussian covariance score等と組み、全方程式について感度H・変動Jと識別性を確認する。既知Cなら係数sampling covarianceは `A^-1 X' W C W X A^-1`、`A=X' W X`。推定CやEIVを含む場合はstacked scoreの `H^-1 J H^-T` を使う。単純な逆Gramへ戻さない。

EIVでは既存の `m=E(X*|Xhat), S=Var(X*|Xhat)` と共有mapping Rを維持する。尤度案Aは `C_beta=C+sum beta_j beta_l R S_jl R'` を使いβ依存共分散も微分する。Bは条件付き平均 `E(y|Xhat)=R m beta` のscoreとして導出できるが、既存lineage loadingに不確かなpredictorを入れる近似が完全な積分と一致するかは追加監査する。EIVに通常REMLを適用しない。説明変数だけから推定した進化rateやsampling covarianceのplug-in不確実性を条件付けるのか、replicate/tipから再推定するのか出力で区別する。

**確定した製品方針**：Bを既定の主解析、Aを共通β・分散・lineage推論の補助解析にする。Bのβ_E検定とAのlineage分散検定を同じ推定量から出たかのように混在させない。Bのモデルベース不確実性計算にAの共分散推定を使う場合も、そのモデル仮定を明示し、Bの推定量に対応するsampling covarianceを計算する。旧 `event` を黙ってAへ置換しない。AとBが異なる効果異質性の例も、Bのβ_Eをtruthとして合格判定する。共通平均モデル下だけの校正でイベント平均推論の対応完了とはしない。

### 3.3 打ち切り観測モデル

既存のobserved lower/upper区間とは別に、全行の**観測機構**を入力する。まず左/右/両側検出限界とexact領域を持つ非情報的censoringを実装する。bootstrapは潜在random effects→未打ち切りY→観測規則の順で生成し、値・ラベル・観測boundsを一緒に再生成する。exactだった行もcensoredになりうる。再fitに古いboundsを渡さない。

interval censoringには、検査時点列・bin境界など全標本空間を覆う規則が必要。一つの観測区間だけから生成機構を推測しない。全区間規則が未指定なら点推定の既存観測尤度は維持できても、unconditional bootstrapは実行不可とする。これは入力不足の対処であり、恒久的な解決策は機構入力と再生成の実装。

固定パターン条件付き解析は代替案にとどめる。正しく行うにはexact値の切断分布だけでなく、共有random effectもパターンに条件付け、条件付きlikelihoodの正規化をそろえる必要がある。censored数やパターンが持つ情報も条件付けで失う。単純なtruncated normal置換は採用しない。

### 3.4 正則化と検定

- `penalty=none`：無罰則marginal likelihood（GLMMはLaplace近似）を最適化し、nuisance-adjusted情報・LR/profileを計算。正則条件が崩れるseparationや境界では通常χ²を自動保証しない。
- 正則化点推定：penalized optimumを残す。無罰則likelihood、penalty、総objectiveを別フィールドにする。罰則曲率逆行列は頻度論的標準誤差と呼ばない。
- **正則化を維持した頻度論的推論の推奨**：帰無制約付きfitからデータを生成し、各replicateで同じ罰則・標準化規則を用いてfull/null双方を再最適化。統計量はpenalized objective差と明記し、経験tail `(1+#T_b>=T_obs)/(B+1)` を使う。これはnuisance plug-inの近似であり、独立校正を要する。
- 信頼集合は候補β0ごとに上の検定を反転する。penalty biasがあるためfull-fit percentile区間や正規曲率区間で代用しない。非連結・非有界集合、root探索失敗はそのまま表現し、便宜的に有限な95%区間にしない。計算予算とMonte Carlo誤差を出力する。
- 無罰則の検定を別fitで行う案も提供可能。ただしpenalized点推定に無罰則の区間を添えて同一推定として扱わず、無罰則MLEと推論を別に報告する。非有限MLEにはpを捏造しない。
- Bayesian案はproper prior、nuisance prior、posterior積分を別途定義する必要がある。現行MAP+Hessianを改名するだけでは採用しない。今回は新MCMC実装まで拡張しない。

lineage heterogeneityも共通null bootstrapへ統合する。半χ²は「単一の識別可能な非負分散・その他のnuisanceは正則な内部・探索なし」等を検証した限定ケースだけにする。joint fixed/variance null、複数境界、少数lineageを自動的に半χ²にしない。bootstrapでは0分散のnullを正確に生成できるよう、微小正数の下限との混同を除く。

## 4. 変更範囲と共通API

共通契約は `estimand / objective_kind / generative_parameters / observation_mechanism / conditioning / statistic / null_constraint / refit_policy` を持たせる。内部名は実装時に型を確定する。生成用共分散をworking covarianceから推測しない。

| ファイル・領域 | 予定する変更 |
|---|---|
| `nwkit/regress.py` | covariance構築のraw/working分離、正規化尤度と均衡scoreの分岐、係数共分散、lineage null、bootstrap、BLUP/予測分散、sampling fraction、omnibus・LOO診断のbasisを整合 |
| `nwkit/gaussian.py` | 正規化Gaussian設定を標準にする。疑似determinant helperの全consumerを点検し、撤去または検証用legacyに隔離。dense/diagonal/low-rank/sparse演算の共通契約を保つ |
| `nwkit/measurement_error.py` | 疑似目的関数の重複修正、β依存共分散のscore/Hessian、共有predictor uncertaintyの生成・条件付け、row scaling除去 |
| `nwkit/phylogenetic_glmm.py` | penalized/unpenalized objective分離、推論適合性検査、null refitと検定反転、観測値+区間の生成bundle、replicated observationへの同機構伝播 |
| `nwkit/regression_pipeline.py` | raw tip→contrast→fitの生成契約、shape再推定、古いboundsの再利用除去、入力マッピングと出力schema。selected contrastの置換で得るtip分布が意図したjoint modelか独立に照合 |
| 新規小モジュール案 | `regression_inference.py` に統計量/制約/実行記録、`regression_observation.py` に観測機構。既存の汎用Gaussian inferenceへGLMM固有契約を押し込まない |
| `nwkit/cli.py`, `regress.py` のCLI転送部、`ordinary_regression.py` | 機構入力・推論mode・デフォルト・全family dispatchを更新。原本のselection callerも追跡 |
| 文書 | `PHYLOGENETIC_REGRESSION.md`, `RECONCILED_SPECIATION_CONTRAST_MATH.md`, `CLI_TSV_CONVENTIONS.md`, `CHANGELOG.md` と専用validation文書。READMEは必要な短いリンクのみ |
| テスト | `test_regress.py`, `test_measurement_error.py`, `test_gaussian.py`, `test_ordinary_regression.py`, `test_cli_contracts.py`, `test_numerical_invariance.py`。新規独立参照・観測機構・null bootstrapテストと小さいCLI fixture |

CLI案：既存 `--event-weighting event|contrast` の意味変更を明示するため `--regression-estimand common|event-average` を導入し、既定を `event-average` にする。既存 `event` はイベント平均の意図を継承し、`contrast` は共通β経路への移行を文書化する。新旧指定の矛盾はエラーとし、raw/precomputed両入口で同じ既定にする。補助解析は明示指定で実行し、主結果を上書きしない。機構入力は `--response-observation-model` と全行のdetection lower/upper column指定を既存のobserved censor boundsから分離する。一般intervalは別のbin/schedule入力へ拡張する。新 `--inference null-bootstrap` と `bootstrap-test-inversion` を区別し、既存 `parametric-bootstrap` の中心化pを新検定と同名のまま残さない。正則化+通常LR/profile/Waldの無条件な頻度論的pを認めない。

結果TSVには最低限 `estimand`, `objective_kind`, `covariance_basis`, `p_value_method`, `interval_method`, `observation_model`, `conditioning`, `bootstrap_attempted/succeeded/failed`, `mc_standard_error` を追加する。行数とevent数を別々に保ち、一般モデルへ機械的に `df=m-p` を付けない。log likelihood/penalty/objectiveは別に出す。旧schemaを読むfigure・集計・examples・CLI schemaテスト・selection callerを一括更新する。

bootstrap失敗を成功例がB個になるまで補充して隠さない。事前にattempt数を定め、失敗理由と成功率を返す。固定B中F失敗ならtailの可能範囲も下限 `(1+exceed)/(B+1)`、上限 `(1+exceed+F)/(B+1)` で記録し、判断に影響する場合は確定pを返さない。seed・入力hash・生成parameter・refitしたnuisanceと固定したnuisanceを保存する。

## 5. 実装後の合格基準（今回は未実施）

### 決定論的・独立参照

1. 生産helperを呼ばないdense NumPy/SciPy式で正規化log likelihood、GLS、REML、EIV covarianceを比較。良条件小行列のobjective絶対誤差 `1e-8`、係数相対誤差 `1e-6` を基本とし、条件数が高い例は先に誤差予算を決める。解析勾配と数値差分、β依存covarianceの微分も検査。
2. 上の期待目的関数反例は修正尤度で真値scoreを `1e-7` 以内、期待最適値を `1e-5` 以内で確認。ML rate比0.95、対応する正規REML比1を閉形式と照合する。Bを採る場合はその不偏scoreとsampling covarianceを別に照合し、MLと同じ期待値を要求しない。
3. 同一行のファイル上のコピーと独立paralog追加を別fixtureにする。Bのloss正規化、単位・factor coding変換、行順序、dense/sparseの一致を検査。Aに独立データ追加時の完全コピー不変性を要求しない。真の観測複製は同じ観測IDまたは正しい共分散で表現する。
4. censoringのCDF確率、exact領域、interval分割の総確率1、zero expected scoreを独立積分で確認。shared latent1次元のGLMMは高精度quadratureでLaplace誤差を測定する。打ち切り0%、片側50%、強い打ち切り、両側、replicate別閾値を含める。
5. Gaussian penaltyの正確なχ²/4例と、別scale/非零βの区間反転を参照にする。Student-tの低次元モデルは直接積分/数値最適化で確認し、Gaussian例だけで合格にしない。
6. 係数・lineage双方で、null constraint、共通罰則、観測機構、測定誤差、shape再fitがreplicateに実際に伝播することをfixtureで確認する。2回bootstrapテストはdispatch/reproducibilityへ改名し、calibration証拠には使わない。

### 独立な統計的validation

開発pilotと最終validationのseed・生成実装を分け、調整後は未使用のseedで判定する。まず小規模pilotで失敗条件を発見し、以下の本検証は実装後に別途実行する。

- event数10/30/100、k=1/2/偏り大、独立eventとcrossed lineage、event/lineage分散0と内部値、balanced/pectinate樹、既知/推定shape、測定誤差0/共有誤差/EIVを代表的に組み合わせる。全直積ではなく主要交互作用を事前登録する。
- censor率0/20/50/80%、左/右/両側/interval、latent variance0/正、既知/未知SD、replicate数と検出限界を変える。80%など非識別条件は「区間を返さない」率も評価する。
- GLMMは少なくともbinomial（非separated/near/complete separation）、Poisson（通常/疎）、censored Gaussian、Student-t/Gaussian/none penaltyを含む。他familyは対応を宣言する前にfamilyごとの観測生成と帰無校正を追加する。
- 主要cellは独立データセット最低2,000。bootstrap検定は少なくともB=1,999を初期値とし、p≈0.05のMC標準誤差約0.0049を記録。検定反転は共通乱数・境界付近の追加replicateでMC誤差による判定反転を抑える。大きな計算費用はpilotから見積もり、実験規模を無断で削って合格扱いにしない。
- 正則・識別可能な主条件の5%検定は棄却率の95%二項区間が事前に決めた許容帯 `[0.03,0.07]` に収まること、95%区間は被覆率の95%二項区間が `[0.93,0.97]` に収まることを基準案とする。境界で離散性により保守的になる条件は上側誤差制御を基準に分け、検出力・幅を併記する。pilot後・本検証前にcell数と多重な合否判定の運用を固定する。
- 成功率99%以上を通常の対応条件の基準案とし、全データに対する返却率、返却区間に条件付けた被覆率、区間なしを非被覆と数えた率を全て示す。選別した成功fitだけの95%被覆を合格にしない。非零βでbias・RMSE、β_Eと共通βの異なるtruth、分散成分bias、lineage検出力も記録する。
- 複数termを探索して報告する場合、上記は事前指定termのpointwise校正と区別する。既存多重検定補正との接続もテストし、選択後のp保証へ拡張しない。

統計基準を満たさないcellは原因をLaplace/識別性/生成機構/最適化/推論に分ける。注意書きだけで対応範囲へ戻さず、正確積分・より妥当な推論へ修正するか、対応条件を具体的に制限する。

実装時には対象pytestから始め、`python tools/check.py full` の全required checksまで実施する。CLI/schema/配布対象変更があればdist検証も実施する。外部依存やskipは明示する。実装後のpytest/full/distの実行状況は実装記録へ分離した。大規模帰無・被覆実験は未実施。

## 6. 他タスクとの境界と実装順序

本担当は回帰の目的関数、EIV、bootstrap、censoring、penalized係数推論、lineage推論を所有する。OUシフトのpBIC・候補探索・選択校正、THRESHOLDのR-hat/ESS、RADTEの区間は変更しない。CLI・文書・共通Gaussian関数を他担当も触る場合は差分単位で統合する。OUのshape固定/推定は回帰内の既存経路まで扱い、shift選択の問題へ拡張しない。RADTEやTHRESHOLDのシミュレーション結果で回帰を校正した扱いにはしない。

原本のGLMM mode polishingは別作業の変更として先に状態を確認し、取り消さない。selectionは予測・CVが主目的なので頻度論的pを追加しないが、fit objectiveや出力のconsumerとして回帰変更の影響を検査する。共有プロトコルを他タスクへ強制するための先行refactorはしない。

1. 最新状態のhash/diffを取り直し、確定済みの「Bが既定、Aが補助解析」に沿ってCLI移行、観測機構入力の仕様を具体化する。小さい独立参照fixtureを追加する。
2. rawな生成モデルCと推論契約を実装し、AのGaussian ML/REMLとEIVを完成。旧determinantとrate流用を除去し、係数共分散・random effect出力をそろえる。
3. Bのevent-average score・stacked uncertaintyを実装し、Aとの推定対象差とコピー不変性を検証する。未完成のBを旧疑似尤度へfallbackしない。
4. 共通null生成/refit基盤を導入し、reconciled係数・lineageとshape再fitをつなぐ。失敗記録・MC精度を実装する。
5. censoring機構入力と観測bundle再生成を実装し、GLMM全refit pathへ渡す。既存bounds-onlyの点推定とbootstrapの入力要件を区別する。
6. 正則化objectiveを分離し、無罰則推論とpenalized null-bootstrap・検定反転を実装。Student-t・separation・Laplace参照で検証する。
7. consumer/CLI/文書/例を更新し、pilot→独立本validation→required checksの順で合格を判定する。

主推定対象は **イベント平均β_Eを既定、共通βを補助解析** とする方針でユーザー確認済み。再確認は不要。残る仕様は「まず全行の検出限界方式、intervalは規則指定時のみ」「正則化点推定は維持し、頻度論的推論は帰無bootstrap/反転へ明示分離」を計画上の推奨とする。この推定対象の合意をもって、大規模validationの実施済み・校正済みとは扱わない。

## 7. 理論参照と今回の実行記録

composite likelihoodでは正規化された周辺/条件付き密度のscoreを組み、感度とscore分散を区別する。通常LRと同一の帰無分布とは限らない。この点は [Varin, Reid & Firth (2011), An Overview of Composite Likelihood Methods](https://utstat.utoronto.ca/reid/research/varin_reid_firth.pdf) の定義・Godambe情報・比統計量の節を照合した。現行目的関数の問題は本文の独立期待値計算に基づく。

今回は原本 `reproduce.py` をimportし、担当4関数だけを実行した。`main()` はresults.jsonを書き換えるため呼ばなかった。worktreeで4反例、原本でcensoring/penaltyの2反例を確認。追加の解析値・一次元数値積分をstdoutへ出し、元results.jsonを変更していない。Python3.10環境でrequestsの依存warningが出たが、4反例は全て計算完了した。
