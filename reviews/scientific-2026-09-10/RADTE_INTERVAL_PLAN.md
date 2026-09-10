# RADTE区間推定：実装・独立検証計画

2026-09-10。対象は科学的妥当性レビュー第7項。今回は読取調査、保存証拠の再集計、既存の小規模テストまで。本番コード変更、commit/push、大規模実験は実施していない。

## 参照状態と確認結果

割当worktree `/Users/kf/.codex/worktrees/573c/nwkit` は開始時にclean、detached HEAD `8d3a4a487bc57b1363c957e1e2bc1e2bb7fab5d6`。ブランチ変更はしていない。AGENTS.mdを確認した。レビュー原本はworktreeにないため、`/Users/kf/repos/nwkit/reviews/scientific-2026-09-10/{REVIEW.md,reproduce.py,results.json}` を読取専用で参照した。原本リポジトリには多数の未コミット変更があるが編集していない。原本の読取時HEADは割当先と同じ。レビュー記載の開始/終了HEADとは異なる。

現在の `radte_studentized.py`, `radte_model.py`, `radte_sequence.py`, `radte_marginal.py`, `radte_uncertainty.py`, `tools/validate_radte_intervals.py` は原本とバイト一致した。ただし保存実験 `examples/radte/interval-coverage-summary.json` のSHA-256との比較では、前3者と検証runnerは一致し、以下は不一致だった。

| ファイル | 現在のSHA-256 |
| --- | --- |
| nwkit/radte_marginal.py | cbbdfe036c3d0596ba8cc02494f61946584dc12ef62f1ccb61746cd93f6a5414 |
| nwkit/radte_uncertainty.py | 9cc8392e49e6e921cee83c3a4507fc34ae921965a5294c58efc4e240dc47ab40 |

したがって以下は保存版の成績であり、現行版のend-to-end被覆検証ではない。レビューのreproduce.py/results.jsonにRADTE計算は含まれず、RADTEの根拠は別のfamily-level CSVである。

`examples/radte/interval-coverage.csv` をPython標準ライブラリで再集計し、本文の分子・分母を確認した。

| 条件 | 方法 | 全family N | 返却 R | 包含 C | C/R | C/N |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 新規4 tips / 2000 sites / SD .3 | Laplace | 200 | 197 | 162 | 82.2% | 81.0% |
| 同上 | studentized | 200 | 199 | 195 | 98.0% | 97.5% |
| 4 tips / 2000 sites / SD .1 | studentized | 200 | 155 | 155 | 100% | 77.5% |
| 4 tips / 2000 sites / SD .6 | Laplace | 200 | 137 | 103 | 75.2% | 51.5% |
| 同上 | studentized | 200 | 199 | 191 | 96.0% | 95.5% |
| 4 tips / 10000 sites / SD .3 | Laplace | 100 | 100 | 79 | 79.0% | 79.0% |

低SDの利用不能45例は全てmarginal推定で、保存された推定SDは全例 `0.0001234098040866 = exp(-9)`。現行marginalのlog SD下限も−9であり、数値下限への到達と整合する。高SDのLaplace利用不能63例のうち62例はGaussian区間の境界越え、1例はactive-bound。低SDと高SDでは対処が異なる。保存CSVだけでは全例のKKT条件や数値下限を外した真の最適解は確定できない。

既存検証は原群200、新規900の独立family、既知SDの再利用200を含む。内部重複100例、16 tips、長配列、SD変化、既知SD対照、複数差分幅によるHessian確認は既に実施済み。これらを新規検証として数え直さない。`python -m pytest tests/test_radte_studentized.py -q` は9 passed（requests依存関係warning 1件）。新規配列推論、MCMCTree実行、全repositoryテストは実施していない。

## 統計的定義と未確定事項

推定対象は固定した遺伝子樹・reconciliation・共有種イベント構造・正しいhard calibrationの下での、指定イベントの年齢。繰返し標本は独立なbranch log ratesと、その枝長から発生した配列で定義する。rhoはまず既知・固定とし、推定rhoの不確実性を暗黙に含めない。置換パラメータ固定の条件付き解析と、毎回再推定する実務pipelineは別に記録する。入力ensembleの広がりはこの頻度論的信頼区間と同一視しない。

- `--uncertainty none` を維持する。Laplaceは無補正曲率、studentizedは局所補正という既存の意味を変えない。
- 各年齢の95%区間は周辺区間。複数年齢の直積を95%同時信頼領域と呼ばない。固定calibrationの一点区間を推論成功率の分母に混ぜない。
- 全試行数Nには入力生成後のfit失敗・区間失敗を含める。返却率R/N、返却条件付き被覆C/R、正しい区間を返せた割合C/Nを別々に出す。最後の量を欠測区間の「被覆率」とだけ呼ばない。
- hard calibration自体が誤って真値を除外した条件では、正しいモデル下95%被覆の合格対象にならない。感度・誤指定実験として、真値除外率と結果を明示する。

branch-onlyでは `y=log(b)-log(d(a))`、既知の正定値Kに対するGLSを局所線形化すると、残差自由度は `n-rank([1,J])`。Kによる白色化はrankを変えないので、rho>0という理由だけで枝数を便宜的な「有効標本数」に置換する案は採らない。ただし非線形性と境界により局所t近似は厳密なpivotではない。

sequence marginalでは `W=V_sequence + sigma² K_conditional` に加えroot pairの非線形積分がある。W全体が未知の単一分散の定数倍ではなく、`n/df` を年齢共分散全体に掛ける操作は配列測定誤差にも作用する。joint-MAPではnon-root枝から推定したSDが曲率計算内で固定される。いずれも線形Gaussianのt公式がそのまま厳密には成立しない。既存studentizedはこの近似を明記しており、現時点で一般的な実装バグとは断定しない。回帰のt区間の前提は [NIST](https://itl.nist.gov/div898/handbook/pmd/section1/pmd131.htm) を参照。

追加で確認すべき点は、(1)低SD境界の統計的非識別と数値制約の寄与、(2)全free ageのうち1座標の不良で全区間を拒否する影響、(3)marginalとjoint-MAPへのauto選択が被覆と欠測選択に与える影響、(4)rate誤差と配列誤差の混合によるstudentizationの過剰/不足補正。条件付き成績を見て事後的にauto選択ルールを変更しない。

## 根本案と代替案

**主案は、未知rate varianceを含む生成モデルに整合した制約付き検定の反転を研究・実装する。** まずbranch-onlyの小次元参照を完成させ、次にsequence marginalへ進む。候補年齢a0で残りの年齢、平均、SDを制約付きで再推定し、LR統計量の帰無分布をrates→枝長→配列という二段階生成と同じ再fit規則で校正する。境界 `sigma=0` を扱える尤度・最適化を用意し、単にlog SD下限を下げたり逆Hessianを正則化して区間を捏造しない。

制約下のnuisance plug-in bootstrapも小標本で保証されない。小次元ではnuisance gridの不利な値でのtail probabilityを参照し、plug-inと比較する。必要ならnuisance confidence set上の上限P値とその集合の非包含確率を組み合わせる保守的方法を検討する。どちらも独立検証で不合格なら正式な推奨区間にはしない。検定反転集合が非連結なら、内側の最初のcrossingだけで閉じない。区間として包絡を出す場合は集合との差と保守性を記録する。

| 案 | 長所 | 制約・採否 |
| --- | --- | --- |
| 現studentizedを維持 | 高速、既存改善の証拠あり | 比較基準として維持。乗数の再調整だけでは境界・配列誤差問題を解決しない |
| rate varianceを含む校正profile反転 | 同じ生成モデルで少数枝・境界・SD不確実性を検証可能 | 主案。費用が大きいため小次元から。joint-MAPの罰則差を通常のLRと呼ばない |
| 局所GLS/REML的な分散補正 | branch-only線形極限で根拠が明確 | sequence混合分散への拡張は導出を要する。自由度や分散成分別補正を独立検証後に採否判断 |
| Bayesian marginal化 | 年齢・SDを同時に積分できる | prior/境界/事後診断とcredible intervalの別契約が必要。既存CIを無言で置換しない |
| 現profile / 現bootstrapへ誘導 | 既存経路を利用可能 | そのまま解決策にしない。profileはχ²閾値、配列bootstrapはsite再標本化で独立rate realizationを生成しない |

境界でのLRは通常のχ²にも一律の50:50混合にも自動的にはならない。対象モデルで導出または帰無校正する（[Self & Liang, 1987](https://pages.stat.wisc.edu/~larget/Stat998/Fall2015/Self-Liang-1987.pdf)）。現配列bootstrapは配列観測変動の感度評価として残すが、全familyにわたるrate変動の校正とのラベルを分ける。二段階生成を新設しても、percentile bootstrapだけで境界解決を宣言しない。

## 不足する検証の優先順位

全因子の巨大直積は組まない。小さいpilotで計算可能性を確認し、主要交互作用を事前指定した検証セルにする。

| 優先 | 条件 | 解決する不確実性 |
| --- | --- | --- |
| P0 | branch-onlyの厳密な小例、既知/未知SD、SD=0と0近傍、正しいhard bound上/直近/内部 | 数値下限、t近似、片側推論を分離 |
| P0 | 同じrate realizationで配列長250/2000/10000/branch-only極限、既知SD対照 | 配列誤差とrate variance推定誤差の分離。対応比較を独立family増加と数えない |
| P0 | rho=0/.5/.9を生成・fit一致、次に生成rho≠fit rho | 相関モデル下の校正と誤指定を分離。root pair・quadratureの精度も確認 |
| P1 | balanced/pectinate/不均衡、浅い・深い内部重複、複数重複、loss、短い内部枝 | 既存の単一nestedケースを超えた自由度・共有年齢・弱識別性 |
| P1 | 固定/幅ありcalibration、境界距離、独立した絶対時間anchorの有無 | calibrationの制約と情報不足。anchorなしは無理に有限区間を要求しない |
| P1 | exact joint-MAP、quadratic marginal、autoを同じ入力で比較 | 自動選択による欠測・被覆選択。外部seedを変えても各fit規則は固定 |
| P2 | HKY/GTR+gammaやcodonで生成し一致/JC69でfit、短配列、欠測・飽和、copy-wide rate shift・重い裾 | 置換/clock誤指定の感度。全モデルの網羅は求めない |
| P2 | 誤topology/reconciliation、外部年代ensemble | 条件付き対象を超える感度。存在しなくなったイベントを成功例から除外せず対応不能率を報告 |

既存の内部重複100例、長配列100例、既知SD200例等は回帰確認用に限定して再利用する。新候補の開発に使ったセル・seedは全てdevelopment扱いに格下げし、独立検証は新しいfamily seed集合、別生成実装、凍結コード/設定で行う。

## 合格基準

1. **独立計算**：二群に分けられるroot重複のbranch-only例では、`theta=log(root duration)` をGaussian二群の平均差として表し、pooled残差とt分布から独立な厳密参照区間を導出する。logit-delta区間との一致は要求せず、その近似誤差を測る。rho既知ではdense KのGLSを参照し、主実装のprecision/gradient関数を流用しない。小配列は全内部状態列挙の尤度と低次元rate積分で検証する。
2. **数値精度**：尺度正規化した小例でlog likelihood差1e-6以内、参照区間端点差1e-3 age-scale以内を暫定基準として事前固定する。積分精度・grid・差分幅を強めても結論が変わらないこと。許容値変更は独立検証を見る前に根拠と共に記録する。rank、root-split不変性、境界KKT、時間単位変更を確認する。
3. **帰無校正**：候補年齢が真値のLR反転試験で棄却率を測る。SD=0/近傍は独立セルとし、optimizer/積分失敗を都合よく落とさない。有限bootstrap回数によるtail確率のMonte Carlo誤差、seed感度、境界付近の数値誤差を報告する。
4. **95%の実用基準**：事前指定した識別可能・正しいモデルの主要セルで、R/Nの片側95%下限≥.95、C/Rの片側95%下限≥.93を暫定採用基準とする。これは「厳密95%」の証明ではなく2ポイントの許容不足を持つ基準。上方の被覆過剰も幅と共に示す。複数主要セルには同時性を調整した判定を使い、単に95%が広い二項区間に入っただけでは合格にしない。
5. **検出力/情報性**：C/Nとその二項区間を必ず併記し、full-domain区間率、calibration幅に対する区間幅、中央値/上位分位、片側率も報告する。常に全domainを返す方法を有用な解決としない。真値と実務上区別したい年齢差を事前に定め、検定反転の排除能力を対応比較する。許容幅は年齢尺度に依存するためpilot後、独立検証前に固定する。
6. **標本数**：pilotは各代表セル20–50 family、最終主要セルは1000程度を起点に精度設計する（95%被覆の標準誤差は1000例で約0.7ポイント）。多重性と返却率を含む必要数を事前計算し、良い結果まで追加する停止規則は使わない。現タスクでは実行しない。
7. **正解との区別**：MCMCTreeとはprior・soft/hard calibration・推定対象を揃えられる範囲だけ比較する。差が小さいことを正確性、差が大きいことを誤りとしない。正確性の根拠は既知生成真値と独立数値参照。MCMCTree側にも独立chainとMonte Carlo誤差が必要。

## 変更範囲と実装順序

1. **証拠と実験仕様を先に固定**：現行版と保存版の差を限定比較し、既存入力の代表例だけ再現。`tools/validate_radte_intervals.py` に実験manifest、生成/fit別rho・モデル・seed、estimator固定、calibration入力、複数target、失敗理由別集計を追加する。現runnerはJC69/1 category/max-age=100固定、target名D一つで、bounds.tsvを渡さない。既存 `simulate` にwidth等があっても現interval runnerでは検証できない点を解消する。
2. **独立参照と原因分解**：新規 `tools/radte_interval_reference.py` と独立生成器を追加する候補。root二群GLS、境界profile、配列誤差の対照を完成させる。低SD45例はactive constraintの種類、log SDのscore、KKT、profile形状を保存する。利用不能率低減は原因ごとに判断する。
3. **生成モデルと境界対応**：`nwkit/radte_model.py`, `radte_marginal.py` のsigma=0極限・制約最適化・識別性を検証し、必要部分だけ修正。`radte_sequence.py`, `radte_sequence_fit.py` に二段階生成/再fitの契約を設計し、root pairの観測可能性を維持する。既存のsite-bootstrapは別方式として保全する。
4. **新区間候補**：`radte_uncertainty.py` で校正profile反転を実装。`radte_studentized.py` は比較基準を維持し、導出に支えられた改善だけを別名または明示版として追加する。小規模検証後に方式を凍結し、大規模独立実験を別工程で実行する。
5. **CLI/出力/消費側を一括更新**：候補CLI `--uncertainty calibrated-profile`（採用前に名称確定）、生成方式・反復数・MC精度を明示する引数を `radte_cli.py` に追加する。`radte.py` のnodes/species/events TSVとmanifestにmethod、nominal level、target別status、境界種別、校正方式、成功/失敗反復数を整合して保存する。target別返却を導入する場合はFitの全体statusとの集約規則を決める。`radte_ensemble.py`, `radte_species.py`, `radte_compare*.py`、図・manifest・TSV利用箇所を検索し、NAと部分返却を一括対応する。
6. **文書・テスト・合否報告**：`RADTE.md`, `RADTE_MATH.md`, `RADTE_VALIDATION.md`, `RADTE_SPECIES_UNCERTAINTY.md`, `examples/radte/README.md` と必要なCLI/TSV規約を更新する。保存証拠は新規ディレクトリに置き旧CSVを上書きしない。`tests/test_radte_studentized.py`, `test_radte_profile.py`, `test_radte_marginal.py`, `test_radte_sequence.py`, `test_radte.py`、ensemble/species/CLI契約テストを更新し、独立参照・境界・失敗時分母を追加する。小規模CIテストと大規模統計検証は分離する。実装時AGENTS.mdとrepository必須チェックを再確認して実行する。

他の3タスク（回帰推論、OUシフト、THRESHOLD）への数理依存はない。共通CLI、CHANGELOG、README、TSV規約を変更する段階だけ担当間で編集を調整する。THRESHOLDのMCMC診断改修をRADTE正確性の代替にしない。species ensemble全体の推論再設計は対象外とし、区間出力の消費契約のみ本担当で整合させる。

実装開始時の推奨判断は「既定none維持、studentizedを比較基準として保存、境界対応を含む校正profileを小次元から開発」。最終採用判断は独立検証後とする。Bayesian方式への全面置換、全条件での有限区間保証、既定方式変更はこの計画には含めない。
