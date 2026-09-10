# NWKIT 科学的妥当性レビュー — 2026-09-10

**結論：科学的な推論として、そのまま信用すべきでない機能がある。** 特に、reconciled 回帰のイベント重み付き疑似尤度、打ち切り応答のブートストラップ、THRESHOLD の収束診断には、今回のレビューで具体的な反例が得られた。OU シフト検出と RADTE には、すでに保存済みの検証結果が示す統計的な弱点もある。

一方、通常の BM/Gaussian 推論、基本的な PIC、PCA、DTT、K/lambda の点推定まで一括して否定する根拠はない。モデルが数学的に定義されていること、実装がそのモデルを計算できること、実データで推定・検定が適切に校正されていることは別々に評価した。

## 対象と証拠の強さ

- 対象は `/Users/kf/repos/nwkit` のレビュー時点の作業ツリー。開始時 HEAD は `cf1e87f31f058d3e40f0b671c7853a2f7416250d`。多数の未コミット変更、新規コマンドも含む。公開リリース全体への断定ではない。
- レビュー中に別作業が進み、終了時 HEAD は DTT を追加する `b8bb512545602fda3e87c29c4c051bf69eaf9102` になった。主要指摘に対応する6ソースの SHA-256 は再確認して一致した。終了時に新たに現れた `nwkit/branch_gaussian.py` と、その後の変更は本レビューの対象外。
- README の機能一覧を入口に、科学的推論を行う主要群の仕様・数式・中心実装・テストを点検した。すべてのソース行、全パラメータ組合せ、全入力形式を証明・網羅検査したわけではない。
- 根拠を「実装の直接確認」「再現可能な反例」「既存実験の再集計・検証」「一般的なモデルの限界」に分けた。既存テストの通過数は、検定の第一種過誤率や区間被覆率の証明とは扱わない。
- 本番コードや既存文書は変更していない。このレビュー、診断スクリプト、結果 JSON のみを追加した。

| 優先度 | 機能 | 判定 | 証拠 |
|---|---|---|---|
| P1 | reconciled 回帰の既定イベント重み付け | 疑似目的関数と推論・生成分布の整合性に問題 | 期待目的関数の反例、コード |
| P1 | `shift` の既定 pBIC | 既知の不正なバックエンドを受け入れる | バージョン検査と保存済み修正検証 |
| P1 | 打ち切り Gaussian の bootstrap | 打ち切りの生成過程を再現していない | 生成器と尤度の反例 |
| P2 | THRESHOLD の R-hat/ESS | 非収束を見逃し、ESS を過大に見せうる | ドリフト・定数列の反例 |
| P2 | 正則化 GLMM の尤度比・profile 推論 | 罰則付き目的関数に通常の χ² 基準を適用 | 数式・汎用推論関数の反例 |
| 科学的制限・強 | OU shift/convergence 選択 | 修正版でも小標本の誤選択が非常に多い | 520 データセットの保存済み証拠を再検証 |
| 科学的制限・強 | RADTE の Laplace 区間 | 小さい遺伝子ファミリーで顕著な過小被覆 | 保存済み family-level CSV の再集計 |

P1 は主要な結論に使う前に対処すべき事項、P2 は該当機能の推論・診断を改善すべき事項。効果の方向や重大性が全データで同じという意味ではない。

## 1. イベント重み付き reconciled 回帰は、通常の ML/REML 推論として正当化できていない

対象：[疑似 determinant と目的関数](/Users/kf/repos/nwkit/nwkit/regress.py:1687)、[係数共分散](/Users/kf/repos/nwkit/nwkit/regress.py:1717)、[境界混合 χ² 検定](/Users/kf/repos/nwkit/nwkit/regress.py:2050)、[bootstrap 生成](/Users/kf/repos/nwkit/nwkit/regress.py:1878)。

同じ種分化に属するパラログの行数を `k` とすると、行固有の共分散を `k` 倍する一方、尤度の標本数を種イベント数に置き換える。さらに determinant は完全な `log|V|` ではなく、イベント平均の分散、または行の周辺分散の平均から作る。二次形式は `r' V^-1 r` のままである。

文書にも composite objective と明記されているので、隠された仕様ではない。しかし、**その仕様に数学的な疑問が残ることは、明記しても解消されない**。

独立なイベントごとにパラログが2個あり、真の残差分散 `s=1`、共有イベント分散 `t=1`、平均が既知の0である Gaussian モデルを考える。真の共分散は各ブロックで `I + 11'`、実装の working covariance は `2s I + t 11'` になる。このとき実装に対応する1イベントあたりの期待負対数疑似尤度は、定数を除いて

```text
f(s,t) = 1/2 * [log((2s+t)/2) + 1/(2s) + 3/(2s+2t)]
```

である。真値 `(1,1)` で勾配は0にならない。20イベントの再現計算では `(-2.083333, -0.416667)`、期待目的関数の最小値は約 `(1.366026, 1.000000)` にある。これは乱数のばらつきではなく、期待値を直接計算した結果である。平均既知の単純化であり、全 CLI 条件の誤差率を測った実験ではないが、定義された生成モデルの分散をそのまま推定しているという解釈への反例になる。

さらに bootstrap は **重み付け後** の `fit['cholesky']` から観測を生成し、同じ疑似目的関数で再 fit する。10イベント・各2行・1つの切片・ML・共有効果なしの場合、再 fit した rate の期待値は元の fitted rate の **1.9倍**になる。反例で元の rate は0.700393、生成分散は1.400785、bootstrap 再 fit の期待 rate は1.330746。重み付けを生物学的な生成分散として再利用しており、自己整合していない。

この目的関数から通常の `(X'V^-1X)^-1` を使い、lineage heterogeneity に `0.5 χ²₁` の P 値を付けることにも一般的な保証はない。正しい composite likelihood でも、通常の尤度比とは異なる参照分布・情報行列が必要になる。[Varin, Reid & Firth (2011)](https://utstat.utoronto.ca/reid/research/varin_reid_firth.pdf)。今回の目的関数では、その前提となる score の不偏性から確認が必要であり、サンドイッチ標準誤差だけで直るとはいえない。

`test_lineage_joint_parametric_bootstrap_reports_calibrated_p_values` は、実際には bootstrap 2回で P 値が0〜1に入ることなどを検査している。[該当テスト](/Users/kf/repos/nwkit/tests/test_regress.py:1278)。このテスト名と文書の “calibrated” は、型I誤差の検証を意味しない。

**対処**：生物学的生成モデルとイベントへの重み付けを分離し、正規化された階層尤度、または導出の明確な推定方程式を定める。変更後に、イベント数・パラログ数の偏り・lineage 効果・測定誤差を変えた独立シミュレーションで被覆率と帰無分布を確認する。`--event-weighting contrast` はこの独自 determinant 経路を避けるが、イベントの影響度と推定対象も変わるため、万能な置換としては推奨しない。

## 2. `shift` の既定 pBIC は、既知の計算不整合を持つ版でも実行できる

対象：[既定 criterion](/Users/kf/repos/nwkit/nwkit/shift_cli.py:41)、[バックエンドの受け入れ条件](/Users/kf/repos/nwkit/nwkit/shift_backend.py:15)。

既定値は pBIC、条件は `kfl1ou >= 3.0.9` だけである。一方、リポジトリ自身の [SHIFT_PBIC.md](/Users/kf/repos/nwkit/SHIFT_PBIC.md) は、元の3.0.9に係数座標と determinant penalty の不整合があり、同じモデルの自由表現と singleton-group 表現でスコアが違うことを示している。

保存済みの estimated-alpha 例では、同じ尤度に対して **32.77624199 と −36.09058668**。修正版は両方約 −36.0905867 にそろう。修正版も同じ3.0.9を名乗るため、バージョン文字列だけでは判別できない。研究用 validation には挙動検査があるが、通常コマンドの入口にはそれがない。

これは「pBIC は有限標本では近似」という一般論とは別の、モデルの表現により criterion が変わる実装問題である。手元の通常 R 環境には kfl1ou がなく、このレビューでは両版の R 推論を新たに走らせていない。元版の不具合自体は保存済み修正証拠に基づくが、通常経路が版番号だけで受け入れる点は直接確認した。

**対処**：修正を識別できる版・機能検査を導入し、既知の不整合がある pBIC 経路を拒否する。単に既定を BIC に変えても、次項の統計的な問題は解決しない。

## 3. OU shift/convergence は修正版でも、現状では探索的な結果として扱うべき

対象：[独立実験の条件](/Users/kf/repos/nwkit/SHIFT_ALPHA.md)、[保存済み集計](/Users/kf/repos/nwkit/examples/shift/alpha-validation/report-data.json)。

8 tips、balanced tree、真のシフトなし、alpha 下限 `alpha*H=1e-7`、two-stage 選択という保存済みの独立 primary 実験で、少なくとも1つの有効シフトを選んだ割合は以下である。

| Criterion | Fixed root | Random root |
|---|---:|---:|
| 修正 pBIC | 47/50 = **94%** | 34/50 = **68%** |
| BIC | 49/50 = **98%** | 49/50 = **98%** |

これは各 criterion が有意水準5%を保証すると約束しているという意味ではない。しかし、シフトがないデータから頻繁にシフトを返すため、選択枝をそのまま適応・収斂の証拠と読むことは危険である。修正済みバックエンドの結果であり、前項のバグだけでは説明できない。

この実験は上限 `alpha*H=10` などを明示した研究用設定で、通常 CLI の既定設定と完全に同じではない。したがって「通常 CLI の偽陽性率は94%」とはいえない。50反復の率にも不確実性があり、非balanced tree・大標本への一般化はしていない。

`verify_shift_alpha_evidence.py` を実行し、**520データセットの再生成、4,160選択結果、160セル、240比較、259,120候補行**を検証した。保存されたパラメータから再計算した尤度の最大誤差は `1.11e-12` 未満だった。これは保存結果の内部整合性を確認するもので、外部 R の最適化を再実行したものではない。

未知のシフト位置を探索する OU モデルで、通常の情報量規準が過剰に複雑なモデルを選びうる問題は理論・実験研究でも指摘されている。[Ho & Ané (2014)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12285)。l1ou 自体には公開された方法論があるが、それがこの fork・小標本設定・収斂拡張の校正を保証するわけではない。[Khabbazian et al. (2016)](https://doi.org/10.1111/2041-210X.12534)。

**対処**：experimental を維持する。実データに近い樹形・種数・ノイズで null と alternative を独立に検証する。alpha 境界の感度、連続最適化、枝位置の非識別性を別々に調べる。bootstrap support や criterion weight を、起源が真である確率と解釈しない。

## 4. 打ち切り Gaussian の bootstrap は観測過程と一致しない

対象：[生成器](/Users/kf/repos/nwkit/nwkit/phylogenetic_glmm.py:2070)、[打ち切り尤度](/Users/kf/repos/nwkit/nwkit/phylogenetic_glmm.py:1054)。

生成器は全観測を正規乱数で発生させ、元データで上下限が指定されていた行を再び欠損にする。打ち切りラベルと限界値は固定され、再生成された値がその限界を超えたかどうかは使わない。元の exact 行は無制約な正規乱数のままである。

検出限界0で、50行が左打ち切り、50行が exact、平均0・SD1の例では、元の exact 行に対して新しく生成した50値のうち **25値が検出限界未満**になった。それでも exact のまま refit される。元の censored 50行は必ず censored のままである。

同じ設定で、bootstrap refit の期待対数尤度の平均パラメータに関する score は、生成平均0で **−39.894228** になる。本来の生成・推定モデルが整合していれば必要となるゼロ score を満たさない。

「打ち切りパターンに条件付ける」という文書の説明だけでは正当化できない。パターンを再生成するなら全観測に適用される検出機構が必要であり、固定するなら exact/latent/random effect を含む条件付き分布と、条件付き解析の定義が必要になる。現在はそのどちらでもない。

**対処**：打ち切り機構を指定できる入出力契約を作り、観測を生成してから censoring を適用する。実装前は当該 bootstrap を利用不可にするか、妥当な推論と誤認しない出力にする。打ち切り尤度そのものが誤っているという指摘ではない。

## 5. THRESHOLD の収束診断は、出力された祖先状態全体を保証しない

対象：[R-hat/ESS](/Users/kf/repos/nwkit/nwkit/threshold_asr.py:262)、[監視対象](/Users/kf/repos/nwkit/nwkit/threshold_asr.py:397)、[判定閾値](/Users/kf/repos/nwkit/nwkit/threshold_asr.py:417)。

R-hat は分割前の通常の chain 平均・分散のみを用いる。ESS は lag-1 相関を AR(1) の式に入れただけで、一般の MCMC に必要な複数 lag の寄与を見ない。監視する量も root liability と可変 threshold に限られ、各内部ノードの category probability/liability を診断していない。

4 chains × 1,000 draws、全 chain に −2→+2 のドリフトを加えた診断用配列で、現在の関数は `R-hat=0.99958, ESS=1092.35` とし、`ok` の数値条件を満たした。単純に chain を前後に分けるだけで `split R-hat=1.35971` になる。さらに全 chain が同じ定数の配列は `R-hat=1, ESS=4000` と判定される。

これは診断関数への反例であって、今回実行した実際の THRESHOLD fit が非収束だったと主張するものではない。修正には rank-normalized split/folded R-hat、複数 lag を使う bulk/tail ESS、定数列を成功扱いしない処理、出力対象の node/category の監視が適切である。[Vehtari et al. (2021)](https://arxiv.org/abs/1903.08008)。

## 6. 正則化 GLMM の罰則付き差分に通常の χ² 検定を適用している

対象：[既定正則化](/Users/kf/repos/nwkit/nwkit/cli.py:3933)、[尤度比と profile 閾値](/Users/kf/repos/nwkit/nwkit/phylogenetic_glmm.py:1881)。

非Gaussian回帰の既定は Student-t 正則化、scale 2.5 である。`likelihood-ratio` / `profile-likelihood` は罰則を含む目的関数を再最適化し、その差を通常の `χ²₁` と比較する。文書にも仕様は書かれているが、通常の likelihood-ratio の帰無分布がそのまま成立するとは限らない。

単純な `Y~N(β,1)` と precision 3 の Gaussian penalty を汎用推論関数に渡すと、統計量は `Y²/4`。帰無分布は `χ²₁/4` である。`Y=2` のとき、実装は `P=0.31731`、この統計量の正しい帰無 tail は `0.04550` になる。この例は Gaussian penalty による関数の反例であり、既定 Student-t の実際の GLMM が必ずこの大きさ・方向で誤るという主張ではない。

十分な情報があれば固定の弱い penalty の影響は小さくなりうる。しかし、この既定値が対象にする sparse/separated data では、まさに penalty の影響が大きい。Wald 共分散も罰則付き目的関数の曲率を使っており、事前分布を使った条件付き不確実性と頻度論的校正の区別が必要である。

**対処**：無罰則 likelihood に基づく推論、明示的な Bayesian 推論、または penalty を含む検定統計量の帰無 bootstrap を区別する。`--coefficient-penalty none` はこの問題を避けるが、separation や Laplace 近似の問題を自動では解決しない。

## 7. RADTE の nominal 95% Laplace 区間には、実測の過小被覆がある

対象：[既存の被覆率報告](/Users/kf/repos/nwkit/RADTE_VALIDATION.md:192)、[今回再集計した CSV](/Users/kf/repos/nwkit/examples/radte/interval-coverage.csv)、[studentized 実装](/Users/kf/repos/nwkit/nwkit/radte_studentized.py:85)。

保存済みの独立 family-level データを再集計すると、Laplace の95%区間は次の結果だった。

| 条件 | 真値を含む / 区間が返る | 返却区間に条件付けた被覆率 | 全 family 数 |
|---|---:|---:|---:|
| 4 tips、2,000 sites、rate SD 0.3、新規群 | 162/197 | **82.2%** | 200 |
| 4 tips、2,000 sites、rate SD 0.6 | 103/137 | **75.2%** | 200 |
| 4 tips、10,000 sites、rate SD 0.3 | 79/100 | **79.0%** | 100 |

配列が長くても改善しない例があり、branch-rate variance を少数の枝から推定する不確実性が重要になる。新規 primary の studentized は195/199、98.0%まで改善したが、別の低SD条件では200例中45例に区間がない。これは条件付きの保守的な被覆であり、全ての family に妥当な95%区間が返るという意味ではない。

この弱点はすでに文書に正確に書かれている。しかも通常の `--uncertainty` の既定は **none** であり、誤った95%区間を既定で強制しているわけではない。studentized も局所 Gaussian 回帰由来の近似で、樹形・reconciliation・置換モデル・校正年齢の不確実性を全て統合したものではない。

**対処**：experimental を維持し、小ファミリーの Laplace を主要な不確実性評価に用いない。studentized の条件付き改善と利用不能率を併記し、研究に近い条件で検証する。MCMCTree は異なる prior/soft calibration を持つため、出力差だけで一方を正解とは判定しない。

## その他の機能の評価

| 機能群 | 今回の評価と残る条件 |
|---|---|
| 通常の BM/Gaussian ASR・PGLS・基本 PIC | Gaussian conditioning/GLS と対応する基盤、独立な数値参照・不変性テストがある。固定された樹形・共分散・観測誤差の仮定下で使うもの。上記のイベント疑似尤度と区別する。 |
| OU・多変量 OU・regime モデル | 数学的モデルは明確。alpha・root・optima・複数分散の識別性と境界、推定パラメータ数に対する情報量が制約。良い fit は自然選択・因果関係の証明ではない。 |
| Mk/HRM/COVARION・stochastic mapping | CTMC pruning/bridge と整合する構成。隠れ状態の非識別性や固定 fitted rate に条件付けた mapping を区別する。図の遷移時刻は歴史の直接観測ではない。 |
| `JUMP-BM`, `MM-BM`, `MM-OU` | fixed-parameter importance sampling として定義され、ordinary ML の IC 比較から外す処理は妥当。高ESSでも未訪問の重要 history は排除できず、独立 seed とサンプル増加が必要。全パラメータを推定した Bayesian posterior ではない。 |
| 個体単位 MV-BM | 進化共分散と個体共分散を分ける観測 Gaussian 尤度は明確。独立個体・共通 within covariance・欠測の仮定が必要。構造的 rank 検査は実用的な推定精度の保証ではない。 |
| PCA | GLS center と進化共分散の固有分解として標準的。軸と固有値の推定不確実性を祖先PC区間は含まない。相関がある入力、単位、ほぼ重複する固有値に解釈が依存する。 |
| DTT | disparity と BM simulation envelope は既存の定義に対応する。帯を観測曲線の信頼帯と呼ばず、MDI を校正済み P 値としない現行文書は適切。 |
| `signal` | K と lambda の基本定義は既存法に対応する。lambda の χ² 近似・境界、小標本、BH の依存仮定は残る。lambda の探索範囲0〜1は明示されており、より広い範囲の他実装と値が違うだけでバグではない。 |
| `regress-select` | training-only preprocessing、nested group CV、post-selection P値を出さない設計は適切。選択頻度はFDR保証ではない。非Gaussian conditional prediction は inverse-link plug-in であり、周辺平均ではない。 |
| reconcile/rooting | LCA duplication-loss や species-overlap heuristic の結果を真のイベントの確定と読まない。ILS、HGT、gene-tree error は仮定外または外部annotation依存。 |
| 樹木の変換・描画・集合演算 | 基本的にはデータ操作。科学的推論の校正という観点では優先度が低い。樹を加工しても生物学的な正しさが増すわけではない。全描画/IO分岐を今回再検証したわけではない。 |

PCA、DTT、signal の方法上の対応は、それぞれ [phytools PCA](https://search.r-project.org/CRAN/refmans/phytools/html/phyl.pca.html)、[geiger DTT](https://search.r-project.org/CRAN/refmans/geiger/html/dtt.html)、[phytools phylosig](https://search.r-project.org/CRAN/refmans/phytools/html/phylosig.html) の公式資料とも照合した。

非Gaussian GLMM 全般では、Laplace 近似の精度は「種数が多い」「疎行列で計算できる」だけでは保証されず、各潜在効果に関する情報量に依存する。これは今回個別に誤差率を測った欠陥ではなく、追加校正が必要な方法上の制限である。[Ogden (2021)](https://onlinelibrary.wiley.com/doi/10.1002/sta4.380)。

## 実行した検証と未実行事項

環境：Python 3.10、NumPy 1.26.4、SciPy 1.15.2。R には phytools 2.3.0、geiger 2.0.11、Rphylopars 0.3.10 がある。通常の R ライブラリには kfl1ou がない。

主要テストは2群で実施し、合計 **676 passed / 23 skipped / 0 failed**。

```sh
python -m pytest tests/test_regress.py tests/test_global_audit_regressions.py tests/test_asr_review_regressions.py tests/test_numerical_invariance.py tests/test_shift_reference.py tests/test_shift_audit.py tests/test_signal.py tests/test_pca.py tests/test_dtt_audit.py tests/test_individual_asr_audit.py tests/test_stochastic_history_audit.py tests/test_regression_selection.py tests/test_radte_studentized.py -q
# 324 passed, 23 skipped

python -m pytest tests/test_gaussian_inference.py tests/test_gaussian_tree.py tests/test_continuous_asr.py tests/test_asr_averaging.py tests/test_full_ou.py tests/test_latent_gaussian.py tests/test_discrete_asr_models.py tests/test_threshold_asr.py tests/test_asr_bootstrap_intervals.py tests/test_asr_comparison.py tests/test_evolution.py tests/test_reconcile.py tests/test_contrast.py -q
# 352 passed

PYTHONPATH=. python reviews/scientific-2026-09-10/reproduce.py
PYTHONPATH=.:tools python tools/verify_shift_alpha_evidence.py examples/shift/alpha-validation
```

23 skips は `NWKIT_TEST_RSCRIPT` が未設定の kfl1ou integration。環境由来の requests dependency warning と、full-OU の有限差分計算中の RuntimeWarning があったが、テスト失敗はなかった。今回の追加スクリプトは Ruff format/check を通過した。

既存 RADTE CSV の family別再集計も行い、本文の分子・分母を確認した。新たな RADTE 配列解析や MCMCTree 推論、kfl1ou の R fit、全 repository の `full`/配布検証、新しい大規模な型I誤差・被覆率実験は実施していない。

再現用 [スクリプト](/Users/kf/repos/nwkit/reviews/scientific-2026-09-10/reproduce.py) と [結果・対象ソースの SHA-256](/Users/kf/repos/nwkit/reviews/scientific-2026-09-10/results.json) を残した。診断用人工配列・解析的期待値を、実データや end-to-end の校正実験と取り違えないよう、それぞれの scope を JSON に記録した。

## 対応の順序

1. reconciled 回帰の生成モデル・目的関数・bootstrap を同じ統計モデルにそろえる。
2. `shift` で不整合な pBIC backend を受け入れないようにする。
3. 打ち切り bootstrap の観測過程を修正する。
4. THRESHOLD の収束診断を改善し、正則化 GLMM の推論契約を整理する。
5. シフト選択と RADTE について、利用条件に沿った独立の誤選択率・被覆率検証を増やす。

単に警告文を増やすことや、同じ実装を使った自己一致テストを追加することでは、これらの問題への十分な対応にならない。
