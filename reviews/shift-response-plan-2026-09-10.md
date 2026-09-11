# OUシフト：pBICの実装識別と選択校正の対応計画

2026-09-10。レビュー第2・3項を統合した**計画のみ**。本番コード、原本のレビュー、別タスクのファイルは変更していない。Rバックエンドの新規実行、bootstrapの大規模再実行、commit/pushも行っていない。

**結論。** 旧IC経路のpBIC受入れを挙動検査で閉じる作業と、既に別タスクで追加されたnative校正経路の独立検証を進める。新しい校正エンジンを重複実装しない。レビュー時点の「通常既定=pBIC」は現在には当てはまらないが、`--selection ic`の既定pBICには指摘が残る。

**参照した状態と所有境界**

- 割当worktreeと原本 `/Users/kf/repos/nwkit` のHEADはともに `8d3a4a487bc57b1363c957e1e2bc1e2bb7fab5d6`。割当worktreeは調査開始時cleanでshiftソースを持たない。原本には多数の変更・未追跡ファイルがあり、shift群も未追跡。以下は原本の未コミット状態を読取参照した計画で、HEADだけでは再現できない。
- 原本の `reviews/scientific-2026-09-10/REVIEW.md`、`reproduce.py`、`results.json` を確認。reproduce.py/results.jsonにはOU選択の新規R再現結果はなく、OU項の根拠は保存済み証拠。レビュー時点は別HEADであり、現在のshift CLIは変更済み。
- 既存タスク「Plan nwkit移植と改善」`01a085c9-487d-71c0-8855-66e503065d8b` の最近の完了内容を読取確認。nativeの探索全体bootstrap、境界処理、700件の証拠は同タスクが追加済み。変更指示は送っていない。実装着手時にその採用状態と差分を再確認し、同タスクの変更を基礎として引き継ぐ。未コミットの所有物をコピーして別実装を始めない。
- `/Users/kf/repos/kfl1ou` はHEAD `63670031d52adf421b5ee6e23968a85f174da4d9`、pBIC修正とテスト等が未コミット。AGENTS.md、Rソース・独立参照テスト・再現スクリプトを読取確認。バックエンド修正の確定・公開識別子はこの所有者側の依存事項。

主要ソースSHA-256（原本の参照状態）：

| ファイル | SHA-256 |
|---|---|
| nwkit/shift_backend.py | `8946565fda853185329b9e702cbd9d42accb7d14941d81bd5d64efc8eb55613e` |
| nwkit/shift_cli.py | `aed23bff31baa4f939191046574182dfe7d63765e0ecdc6a7ed2e17f625bab50` |
| nwkit/shift_calibration.py | `05651f9975884c9cbdadf360497fc6a1fa09e52b9130e83932d76348bb965fdc` |
| nwkit/shift_candidates.py | `87f2f98cc05deaa848b0d573d95cbcea41c6d466b77a8783a78c3449f36094ef` |
| nwkit/shift_calibrated_output.py | `3937df8036b687c3b57adb11d453580ff07680c25a8adc6577ccb95602add4b9` |

**終了直前の並行変更**：上表は調査中に読んだ状態である。最終確認時に原本のcalibration本体が `e30e9aab8c4886811c315817e0b730e7d3966bfddb743899ec18c6cded6809cf`、候補列挙が `bdf10b084ea6d9b8cda0ffe0a7758b617a41e694ad481aed01ff7c65f4e510a0`、出力が `0bc3797f92f7102fc56772bf62aaa69522571aa036d9a98530b436922c6711ea` に変わった。保存snapshotとの読取diffで、正規方程式からQRへの変更、bootstrapの64件ごとの処理、入力検査追加、convergence-off候補の事前除外を確認した。これは別作業の変更で、本計画ではその結果同値性を再実行検証していない。従って「protocolのソースhash一致」はその変更前に確認した事実であり、保存700件の証拠が最新実装を直接検証するという意味ではない。probe入口とCLI既定の確認対象hashは最終確認でも同じだった。実装開始前に最新差分と新たな証拠の対応を確定すること。

**確認できた問題、既に変わった点、未確定な点**

| 論点 | 今回の判定 | 次の対応 |
|---|---|---|
| 元版/修正版とも3.0.9 | 確認。通常R adapterは現在も最低版だけを検査。研究driverの固定αprobeはadapterで使われていない | pBICを計算する入口すべてで修正能力を確認 |
| 座標不整合 | 行列変換と保存CSVが整合。自由表現の推定αスコア32.77624に対しsingletonは−36.09059、修正後差は1.42e−7 | Rで両インストールを分離して受入・拒否を確認。自由/singletonの一致だけで全pBICを認証しない |
| 旧実験の過剰選択 | 保存証拠の再検証を今回実行して通過。8-tip balanced、αH下限1e−7・上限10、two-stageでpBICは47/50と34/50、BICは両rootで49/50 | 通常CLI既定の率と呼ばない。IC、探索法、境界、最適化を切り分ける |
| 現在の通常CLI | `--selection calibrated`、4–16 tips・最大2シフト、B=199、level=.05、convergenceは既定off。ICは明示選択 | 旧研究文書の「production unchanged」は当時の履歴と明示し現行仕様へリンク |
| 新校正の700件 | 元ソースhash一致と保存行を今回再集計。帰無10/250。別保存の既定convergence-offも10/250 | 改善の証拠だが一様5%保証ではない。旧520件と異なるデータなので94%→4%をpaired効果としない |
| セル別校正 | 新primary固定rootは5/50、16-tip固定rootは3/25。小標本で率が粗く、過大とも妥当とも断定できない | pooled平均だけで合格させず、新seed・新樹形でセル別評価 |
| α・分散推定 | nativeはαHの0/∞と有限25点、既知SEありは分散49正点と0のグリッド。連続最適化ではない | 選択統計量の定義にグリッドを含め、連続参照と感度評価 |
| 収斂 | nativeは共有効果を含む候補を共同列挙。旧ICは枝選択後のbackward merge | 共同探索の完全性と、収斂が識別できるかを別検査 |

今回 `verify_shift_alpha_evidence.py` が確認した範囲は520入力の再生成、4,160選択、160セル、240比較、259,120候補行、保存パラメータでの尤度最大誤差 `1.1032e−12` 未満。R最適化の再実行や大域最適性の証明ではない。新証拠の再集計では、`records.jsonl.gz` SHA-256 `b4ed1b100b9e1a0cc020d177541fdae4ce3a70ad12763014221d1c91fc535f0e` とprotocol内5ソースのhash一致を確認。新700件のfitそのものは再実行していない。

**採用する統計的定義**

1. 有効シフトは最も近い選択祖先と異なる群への変化。帰無誤選択は「少なくとも1有効シフト」。保持した枝の数とは分ける。枝の完全一致、末端群分割のラベル不変一致、末端期待値のRMSEを別々に評価する。
2. pBICは明示した係数座標・位置ペナルティ・推定共分散パラメータ数に従うモデル選択規準であり、検定P値ではない。自由/singletonは同一モデルの再表現なので同じ定義のpenaltyになるべき。任意の再パラメータ化で先験体積を無視して不変になるという主張ではない。収斂拡張のpenaltyの導出は通常pBICと別に点検する。方法論の一次資料は [Khabbazian et al. (2016)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12534)。同論文の検証条件をこのfork・小標本・拡張の保証に読み替えない。
3. nativeは `z=Qy, Q1=0, QQ'=I` のGaussian contrast密度を使用する。フルデータMLやモデルごとに平均を消去するREMLとは区別する。`z ~ N(QXβ, Q(vK+D)Q')`。Dは既知の独立観測誤差、vは過程の末端分散。ultrametricで同じα・σ²ならfixed/random rootの共分散差は共通rank-one項でありQが除去する。今回8-tipの直接計算でその投影残差は `5.63e−18` 未満だった。root不確実性を推定したことにはならない。
4. 現行検定を正確に記述する：入れ子族F0（0シフト）、F1（高々1シフト）、必要ならFshared（高々1非baseline効果）、Fallを用い、`Tj=2(max_Fall ℓ−max_Fj ℓ)`。Fjの選択fitから生成し、各replicateで全候補と全nuisance探索を再実行、`p=(1+#{Tb≥Tobs})/(B+1)`。最初の非棄却族で停止。既知誤差なしの生成scaleはRSS/(d−q)、尤度のscaleはRSS/d。この自由度補正も事前固定した手続きの一部として、ML scaleとの比較を行う。
5. no-shift帰無で何かを選ぶには第1段を棄却する必要がある。しかしplug-in nuisance推定と有限Bがあるので、この事実だけで正確な5%とはならない。後段の族選択は収斂の真偽や各枝の存在に対する個別の5%検定ではない。推定された共有効果の採用を「収斂を証明」としない。
6. α=0は分けて定義する。有限のOU最適値を保ってα→0ならシフト効果は消える。現nativeはscaled effectを固定するドリフト極限まで候補に含めるため、有限θのOU族の単純なBM置換ではない。α=∞では末端平均・vを扱えても有限αとσ²の分離はできない。α境界を支持する場合の有限θ/optimum effectはNAを維持する。`optimum_identifiable=true`も現在は診断規則の通過であり識別性の証明ではない。
7. bootstrap supportは、選択済みfitを生成モデルとした再選択頻度。選択の安定性を表し、枝が真である確率・P値・区間被覆率ではない。criterion weightも事後確率とはしない。群分割支持と枝支持、試行/成功/失敗数と分母を保存する。現段階でモデル平均化を解決策にしない。

**根本修正と代替案の比較**

| 層 | 採用案 | 代替案と不採用・条件 |
|---|---|---|
| バックエンド受入 | 正しい修正版識別子＋公開APIの小さい挙動probe。pBIC使用前に失敗を検出し計算を拒否 | 最低版を上げるだけでは同版改変・異なるlibraryを識別できない。hash固定だけは将来の正しい版まで拒むため再現用lockに限定 |
| pBIC修正 | kfl1ou側の座標変換、実際のconstrained αでの情報行列、固定α等のパラメータ数修正を確定 | Python adapterに補正式を後付けしない。BICへの自動置換も誤選択の解決にならない |
| 小規模樹の選択 | 既存native全探索bootstrapを基礎にし、独立校正の合格範囲を明示 | pBIC/BICを5%検定として使わない。大きい樹への自動fallback・外挿もしない |
| nuisance校正 | まず現在のplug-in法を事前規定した条件で監査 | 不合格ならnuisanceの信頼集合上の最悪tail＋集合外確率補正等を検討。ただし有限グリッド上の最大値だけで連続空間の保証を主張しない。二重bootstrapは計算費用と追加近似を評価してから |
| α最適化 | 研究用の独立連続profileを先に作り、現行グリッドとの差を測る | 下限を恣意的に上げない。連続探索を採用するなら、観測だけでなく全bootstrapにも同じ探索を適用して校正し直す |
| 非識別性 | 平均予測・scaled effect・候補同値性を出力し、有限θと区別 | 警告だけで巨大θを確定推定として返す、同じtip partitionだけで異なるOU履歴を統合する、は不可 |

座標補正は `m_b=w_b δ_b, w_b=1−exp(−α parent_age_b)` に対し `log|Cov(intercept,δ)|=log|Cov(intercept,m)|−2Σlog(w_b)`。今回任意の正定値3×3行列で変換を直接確認し残差は `8.89e−16` 未満。

さらに固定root・非退化の自由シフトdesignでは、小さい正のαでoptimum座標の情報行列に `2k log α` が現れる。固定データで平均変位を再fitした尤度が有限でもpBICは低下し得る。8-tip保存fixtureの独立dense計算ではα `.01, 1e−4, 1e−6` で推定αとして数えたpBIC式は `14.3447, −4.0652, −22.4858`、尤度は約 `−7.50` に留まった。これは**同じ規準をαに沿って評価した数理診断**であり、backendがαをpBIC最小化で選ぶという意味ではない。backendのα推定は尤度fit経路に依存する。α=0でshifted pBICを∞とするガードだけでは小さい正のαの問題を解消しない。近似の非正則領域として扱い、正則な条件を仮定したpBICの較正を別に評価する。

**変更対象と契約**

| ファイル群（原本基準） | 実装時の変更内容 |
|---|---|
| `nwkit/shift_backend.py`、必要なら共有probeモジュール | 検査は実際にfitする同じR executable/library内で実行。固定αの座標・推定パラメータ数、推定α、constrained再fitを区別。NaN/∞/欠落も失敗。probe ID/version、実行path、解決library、package識別/hash、スコア残差を返す |
| `nwkit/shift.py`, `shift_results.py`, `shift_cli.py` | IC pBICの未検証実行を拒否。BIC等を勝手に代用しない。通常設定と研究設定、実際のboundをJSONへ。calibratedとIC限定optionの混用拒否を検査。probe失敗時に成功らしい最終結果を残さない |
| `nwkit/shift_calibration.py`, `shift_candidates.py` | 別タスクの実装を維持して監査。α/variance profile、候補族、rank判定、tie policyを検証。連続探索を採用する判断が出た場合だけ同じエンジンを全再標本に使用 |
| `nwkit/shift_calibrated_output.py` | 族の意味、ML/生成scale、有限grid、fit失敗方針、seed/B、boundary、同値候補と識別診断を再現可能に保存。既存schema 5/6の消費者と対応を確認してからschema変更。非有限値はstatus＋null/NA |
| `nwkit/shift_convergence.py`, `shift_bootstrap.py` | legacy二段階とnative共同探索を区別。bootstrap再fitにも同じbackend検査・設定を適用。支持率の分母、失敗、群分割と枝を分離 |
| `tools/shift_alpha_backend.R`, `validate_shift_alpha.py`, `validate_shift_joint*.py`, `shift_joint_backend.py` | productionと共通の能力契約を使用し、研究driverだけに検査を閉じ込めない。旧二段階/共同のpaired比較は同じ候補範囲・bound・root・error・連続最適化を用いる |
| `tools/shift_simulation_cases.py`, `validate_shift_calibration.py`, `shift_alpha_design.py`等 | 開発データと独立検証を分離し新protocolを先に凍結。通常CLI-default laneと研究laneを別保存。独立生成器と独立dense参照を本体helperに依存させない |
| `verify_shift_alpha_evidence.py`, `verify_shift_calibration.py`等 | 旧証拠は不変。新verifierの読取監査と再fit/出力を分離する。現verify_shift_calibrationは入力directoryにaudit/recordsを書き込むため、そのまま原本に実行しない。生成状態・CLI wrapperのhashも保存対象に広げる |
| `tests/test_shift*.py`、新backend能力テスト | original/fixed R環境、契約、独立数式、境界、探索・支持の意味をテスト。Rなしのmockだけで合格にしない |
| `SHIFT.md`, `SHIFT_PBIC.md`, `SHIFT_ALPHA.md`, `SHIFT_JOINT.md`, `SHIFT_CALIBRATION.md`, `SHIFT_VALIDATION.md` | 現在の既定と歴史的研究条件を明示。有限θとdrift、stagewise P値とsupport、full likelihoodとcontrast densityを区別 |
| `examples/shift`、README/ASR.md、CLI契約テスト、MANIFEST.in、tools/check_dist.py | 新例は別directory。help/例/JSON/TSVを同期。ASRはregime-mapを再fitするので、境界モデルや選択不確実性を受け渡したと扱わない。READMEは短い入口に留める |
| kfl1ouのR/shift_configuration.R, R/convergent_regions.R、対応Rテスト・文書 | 別所有タスクの修正確定に依存。こちらから編集しない。受入probeの仕様と正式識別子を確定後に接続 |

probeは「自由/singleton両方に同じ誤りがある」場合を見逃せるため、通常入口の小fixtureにも独立dense参照の期待値または情報行列を加える。詳細行列計算はテストで網羅する。キャッシュする場合はR pathだけでなく実際のlibrary内容・probe版をキーにし、別インストールや上書き後の誤受入を防ぐ。テスト用隔離libraryは再現性のためで、利用者のlibraryを変更しない。

**合格基準（実験開始前に凍結する提案）**

1. **実装整合性**：元3.0.9の既知不整合fixtureを必ず拒否、修正版を受入。BICは非影響fixtureで変化なし。独立dense尤度/平均は絶対誤差1e−6以下、well-conditioned fixtureのpBIC/情報penaltyは2e−6以下。悪条件では単純な緩和で通さず、condition数・高精度参照・失敗状態を検査。fixed/random root、α固定/推定、nested/disconnected/shared、SEあり/なしを含む。
2. **探索**：8/16-tip研究候補は195/899と独立の列挙で一致。候補ID一意、欠落なし、失敗を除外したまま「完全」としない。jointの改善は同じ連続設定・含有関係が確認できるときだけ比較。離散全列挙を連続大域最適と呼ばない。新nativeについてもconvergence-off/onの族と候補数を別検証する。
3. **境界・数値**：branch単位生成/密行列からfinite α、0、∞、v=0、既知SE不均一、深い非balanced treeを検査。tip並べ替え、定数加算、trait/time単位変更で必要な変換則が成立する。dense連続profileをログgridの各局所極大から精密化し、粗いgridでの選択差を記録。連続モードを導入する場合の小fixtureで尤度差1e−6以下を要求。有限gridの現行モードは近似として表示し、独立検証で連続参照との選択不一致率の95%上限1%以下を暫定許容案とする。不合格ならgridを再設計し新seedで再校正。
4. **帰無校正**：1データセットを独立単位とし、同一データの8設定やbootstrap反復を独立標本数に数えない。通常CLI laneはconvergence off/on、4/8/16 tips、balanced/非balanced/浅い分岐、αH=0を含むweak/中/大、SE=0/既知/不均一を事前指定。random-rootのα=0はstationary分布未定義なのでその生成セルを作らず、contrast/BM極限を別記する。長枝・短枝・v≈0、true nullを含む。研究laneは旧ICのboundと探索法を固定して比較する。
5. **反復数・判定**：主要セルはまず独立1,000 datasetsを計画し、pilotで計算量だけを見積もる（本タスクでは開始しない）。最終Bは999以上を基本に、CLI既定B=199も別に検証する。名目.05に対し、主要セルの誤選択率の同時片側95%上限が.075以下を実用上の暫定合格基準とする。これは5%の厳密保証ではなく許容差2.5ポイントの工学的基準。Bonferroni等でセル数を反映し、必要反復数は設計段階で再計算する。単に区間に.05が入る、全体平均が低い、だけでは合格にしない。失敗は分母に残し、全失敗を誤選択と数えた上限も基準を満たすこと。セルを見てから追加して停止する場合は別の逐次設計を先に定める。
6. **対立・収斂**：single/distinct/shared、nested return、近接した異なる最適値、非識別な履歴を含める。検出力、過小/過大選択、枝位置・末端群分割、有限θ返却率、平均RMSEを別報告。distinctをsharedと断定する率も評価する。校正改善の代償を隠さず、強い識別可能な事前指定セルで現行nativeより検出力を5ポイント超悪化させないことをpaired独立データの区間で確認する。弱いpullで高い回復率を必須とせず、推定不能を正しく返すことを重視する。
7. **区間・支持の解釈**：現α profile cutoff 1.920729は診断であり95%CIと呼ばない。今回の範囲で新しいCIを作る必要はない。将来CIを出す場合は選択過程も含めて校正し、全生成例に対する返却率、返却条件付き被覆率、全例での「返却かつ被覆」、区間幅を別報告。名目95%をうたうなら主要セルの被覆率同時下限.925以上等を先に設定する。bootstrap supportについてはtrue/false枝ごとの支持率分布を調べ、0.95支持を95%真確率と読み替えない。

**実装順序・依存関係・必要な判断**

1. 原本と既存タスクの状態を再確認し、採用する未コミット差分/ソースsnapshotを確定する。今回のhashから動いていた部分は再読する。回帰・THRESHOLD・RADTEの推論実装は担当外。共有cli.py/provenance/文書を変更するならその所有者と範囲を調整し、他のデバッグを待たずshift固有の契約テストを先に進める。
2. kfl1ouの既存修正を基にprobe契約と独立期待値を確定し、NWKIT通常IC入口の拒否とprovenanceを先に実装する。修正版公開が未完でも、元版を拒否する挙動・小テストは独立に進められる。正式release番号はここで捏造しない。
3. 現nativeについて平均族・境界・段階的検定の定義を固定し、独立参照・候補網羅・CLI/ASR出力契約を検証する。α=0のdrift候補を既定のOU探索に含めるか、明示別モードにするかは科学的推定対象の判断事項。推奨は現実装の機能を保ちながらモデル族を明示し、通常有限OUとの結果を別扱いにすること。
4. 小さい連続profileとnuisance感度比較を行い、固定grid継続/適応的連続探索の判断をする。数値精度改善だけを先に観測fitへ入れず、bootstrap生成・再fitへ同時反映する。
5. 独立protocol・新seed・合格基準・計算予算を凍結してから上記帰無/対立実験を実施する。nativeの現行既定を基礎とするが、現700件だけをもって科学的検証完了や5%保証とはしない。非合格セルが残ればその範囲を研究用途とし、nuisance校正改良後に新しい検証群で再評価する。
6. schema・CLI・例・文書を同期し、対象pytest、実R両環境、interface/CLI契約、ruff/mypy、配布検証とリポジトリ所定チェックを実行する。別タスクが報告した通過数を今回の検証結果として転記しない。本番コードの変更はこの計画作成では行わない。

実装へ進む際の主な判断は、**drift極限を含む推定対象の採用、5%に対する実用許容差と検証予算、kfl1ouの正式な修正版識別契約**の3点。BICへの切替、α下限の引上げ、警告追加、自己一致テストのみで完了とはしない。
