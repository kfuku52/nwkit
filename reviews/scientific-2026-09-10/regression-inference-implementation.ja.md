# 回帰推論の実装記録

2026-09-10。レビュー第1・4・6項への変更。割当worktreeの未コミット変更として実装した。

## 実装した契約

- reconciled Gaussianの既定はraw contrastのイベント平均β_E。W_ii=1/k_e、L=(X'WX)^-1X'Wで推定し、共分散は補助モデルCを用いたLCL'とする。EIVではCをβ_Eで評価する。共通βの階層尤度は`--regression-estimand common`で選ぶ。旧event/contrast指定との矛盾を拒否する。
- 生物学的共分散をパラログ数で膨らませず、分散成分を実際のnと完全なlog determinantで推定する。event pseudo-determinantは拒否する。event平均の自由度を機械的にイベント数−pとはしない。
- bootstrapは生物学的Cから生成する。shape再fit時はtip空間の進化・sampling noiseを保持し、共有効果を追加することで、選択外contrastや祖先値とのsampling cross-covarianceを保つ。
- Gaussian warm startが分散境界に張り付く問題を修正し、各データセットで内点からも最適化する。全モデルが帰無モデルより悪い目的関数差をゼロへ丸めない。
- censored Gaussianは全観測の検出限界・interval partition・uncensored機構を明示入力し、潜在応答から値・打ち切りラベル・上下限を再生成する。観測boundsのみのfitは可能だがbootstrapには機構が必要。
- 既定のStudent-t正則化は点推定として維持し、罰則曲率に通常の頻度論SE/p/CIを付けない。`null-bootstrap`では各係数の制約付き帰無から生成し、同じ罰則の全モデル・帰無モデルを各データで再fitする。
- 係数profileは指定グリッド上で各帰無を検定して反転する。p、Monte Carlo SE、採否をJSON列に保持する。未評価点を補間した連続CIを返さない。無罰則profile-likelihoodの端点探索失敗をWald区間へ置換しない。
- bootstrapの試行数を固定し、失敗データを別データで置換しない。失敗がある場合は推論をエラーで終了する。成功結果に試行・成功・失敗数、推定対象、目的関数、尤度帰属、区間・pの方式を記録する。
- 旧cluster-HC1は従来の分散標準化lossを残し、`legacy-variance-weighted-sensitivity`と明示する。lineage/random effect推論も補助共通βモデルに帰属させる。

## 検証

独立参照テストは、2イベントのβ_E=2と共通β=13/11の差、dense/diagonal/low-rank/GMRFのLCL'、既知二次モーメントによる分散成分s=t=1の回復、ML bootstrap rate期待値0.95とREMLの1、censoring期待score=0、正則化Gaussianのχ²/4帰無分布を検証する。4,000回の一次元乱数は小さい解析参照の確認であり、GLMMの4,000データセットfitではない。Student-t_3罰則についても、一次元正規likelihoodの直接最適化と独立な正規分布の裾を2,000回の小さい乱数参照で照合した。これはStudent-t GLMMの厳密積分検証とは異なる。

Poissonの正則化null-bootstrap、censored Gaussianのreplicate-aware CLI→観測生成→再fit、separation時のpenalized null経路について、小さいB=2の統合テストを実施する。これらは配線・再現性テストであり、5%検定や95%被覆の証拠とは扱わない。

- 新規推論契約テスト: 27 passed。
- 最終メタデータ・最適化変更後の回帰系テスト: 115 passed。
- ordinary regression既存スイート: 93 passed。追加のcensored CLIケース: 2 passed。
- lineage最適化修正とCLI契約: 48 passed。
- Ruff lint/format、mypy（164 source files）、pip check、Bandit、pip-audit: 成功。既知脆弱性なし。
- 複雑度チェック: hard limitを全て通過。上限緩和や例外追加なし。
- `python tools/check.py dist`: wheel/sdist内容・再現性検証に成功。
- `python tools/check.py full`: 終了コード0。全体pytestは **3,098 passed / 27 skipped / 6 warnings**（908.95秒）。coverageは **84%**（必須80%を通過）。全体実行開始後に追加・最終調整したケースは上記の個別テストでも確認した。スキップ27件は成功へ数えず、スキップされた経路の実行検証は未実施。

Python環境のrequests依存warningは発生するが、pip checkではbroken requirementなし。全体実行には既存ASR比較のpandas FutureWarning、full OU/RADTEのSciPy数値最適化warningも含まれる。個別スイート件数は重複を含み、合計して独立テスト数とはしない。

## 計画から具体化した範囲と未検証事項

- β_EのLは共分散nuisanceに依存しないため、その一次微分はゼロ。実装はmodel-based plug-inであり、説明変数モデルの推定誤差や任意の平均異質性・モデル外依存性に対する一般的な頑健被覆を保証しない。
- Gaussian係数bootstrapはfull-fit中心化pとpercentile区間の近似を維持し、正しい生成Cと選択された推定量へ修正した。GLMMの係数別帰無検定とは区別する。
- 正則化した連続CIの自動端点探索は実装せず、評価済み候補グリッドに限定して検定反転を返す。グリッド内の間隙・外側・非連結性を有限区間で隠さない。
- 大規模な帰無棄却率・被覆率実験、Student-t GLMMの独立な厳密積分比較、全family/境界/shape/EIVの科学的校正は未実施。計画中の2,000データセット、B=1,999等の合格基準を満たしたとは主張しない。
- Laplace近似、nuisance plug-in、離散bootstrapの精度限界は残る。新しい推論経路の実装済みと、全対応条件で校正済みとは区別する。

## 本体への統合

`master`の9e2dd93（GLMM multistart改善）へ統合。既存のscalar mode polishingと無罰則multistartを保持し、他タスクの未コミット変更はコミット対象から除外した。

本体に新設された校正ツールの旧pseudo-likelihood引数も更新し、補助共通係数Gaussian MLの帰無検定と明記した。保存済み実験のソース・結果は変更せず、旧実装の校正結果を今回の推定量の検証に流用しない旨を文書化した。

統合検証は推論契約・ordinary回帰・CLI・校正ツールの193ケース。初回192成功と旧引数による1失敗を確認し、修正後は校正ツール20ケースが全て成功（重複を含む）。Ruffとmypy（200 source files）も成功。全体チェックとdistの実施結果は上記worktreeでのものを保持し、本体では影響箇所の統合検証を追加した。

終了判断: 指摘への実装修正と本体統合は完了。連続区間の自動端点探索と新実装に対する大規模な頻度論的校正は未完了であり、科学的検証全体が完了したとはしない。このタスクは実装タスクとして終了可能。
