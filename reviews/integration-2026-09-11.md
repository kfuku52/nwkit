# Integration review — 2026-09-11

## 対象

レビュー開始時の本体は `master`、HEAD は `0caec89`、更新後の
`origin/master` は `cf1e87f`。未 push の23コミットを対象とし、開始時の
本体と実在する別 worktree に未コミット変更はなかった。
修正は本体で行い、既存のパッチ更新方式に従って `0.43.16` に更新した。

コード差分、CLI 契約、数値境界、入出力、依存関係、全テスト、配布物を確認した。
全入力・全外部実装について正しさを証明するものではない。

## 検出した問題と修正

| 問題 | 影響と修正 | 検証 |
|---|---|---|
| 定常固有値の負の丸め誤差 | RADTE の長い枝で遷移確率の行和がほぼ0に崩れ、尤度計算が失敗した。定常固有値を数理的な値0に固定した。 | 独立な行列指数・定常分布との比較、全対応置換モデル、尤度勾配、profile LR の再推定。 |
| ASR の同じ定常成分の問題 | 対称 Mk 過程にも適用。分離した状態群は個別に固有分解し、速い群の丸め誤差が遅い群の固有値を上回る場合も保持する。 | 連結・非連結グラフ、速度比 `1e20` の独立 JC 閉形式、ASR/Mk 関連117テスト。 |
| 空の SciPy 境界 | 年代が全て固定された枝長ベース RADTE で、新しい SciPy が空の `Bounds` を拒否した。推定する nuisance parameter がない場合は、その境界を構築しない。 | RADTE CLI、固定年代出力、calibrated profile、通常 profile。 |
| 微小な偽の正の分散 | Python 3.10/SciPy 1.15 の SLSQP が真のゼロ分散を `7.8e-17` 程度として返し、interval の境界診断が変わった。尤度が丸め誤差範囲で同等の場合にゼロ分散を明示的に再推定し、右微分も確認して採用する。 | 実際の RADTE 推定、既知のゼロ／微小な正の最適値を持つ独立な二次目的関数。 |
| pandas の文字列列 | 未計算の Monte Carlo SE を空文字列で初期化すると、pandas 3 で bootstrap の数値を書き込めなかった。数値の欠測値として全該当 producer を統一した。 | 自動進化パラメータを再推定する実際の regression bootstrap。TSV の欠測表示は維持。 |
| root の互換 alias | `--candidates-out` に共通規約の `--candidates_out` がなかった。alias を追加した。 | 全サブコマンドのインターフェース契約。 |
| 一回限りの iterator | `ShiftLayout.build` が検証時に枝 ID の generator を消費し、指定シフトを無言で失った。最初に tuple 化した。 | 枝ごとの独立な平均伝播との一致、iterator 内の重複拒否。 |
| locale 依存の読込 | simulation JSON と predictor 一覧が非 UTF-8 locale で日本語名を読めなかった。UTF-8 を指定した。 | 非 UTF-8 locale を模した実際の CLI、simulation 出力、nested regression selection。 |
| 検証資料の厳密な浮動小数比較 | 異なる NumPy/SciPy/BLAS で `1e-16` 程度の差を資料破損と誤判定した。生成配列の形状は厳密に保持し、計算値だけを丸め誤差範囲で比較する。 | 保存済み240データセット・480 fit の読取専用再実行、source hash、seed、件数、p値、採否の厳密一致、改変拒否。 |

生成配列は `rtol=1e-13, atol=1e-14`、LR 統計量と再計算した同時信頼上限は
`rtol=1e-12, atol=1e-12`。その他の構造・値は従来どおり厳密比較する。
保存済みの観測、結果、source snapshot、summary は変更していない。

## worktree と取り込み漏れ

- `/private/tmp/iqtree3-push-nwkit` は clean、先端 `cf1e87f` は本体の祖先。
- 登録だけ残った `/private/tmp/nwkit-selected-deps.XWihmZ` の先端 `5bb4809`
  も本体に取り込み済み。ディレクトリは存在しない。登録情報は削除していない。
- Codex worktree 保存先と `/private/tmp`・リポジトリ周辺も確認した。
  生存する追加の NWKIT worktree は見つからなかった。
- stash は空。古いローカルブランチの独立した packaging commit は
  `git cherry` で本体とのパッチ同値を確認した。
- 到達不能コミットも確認した。旧 PGLS/ASR 変更はパッチ同値、旧配布物修正は
  バージョン更新込みで統合済み、古い依存更新は現行の更新に置き換えられている。
  過去の convert stash も本体に統合後、追加機能を含む状態に更新されている。
- archive snapshot `ad9108c` と `649b617` の変更ファイルを本体と照合した。
  実装は統合されていたが、後者のレビュー資料4ファイルが本体に欠けていた。
  `SHIFT_RESPONSE_VALIDATION.md` は、そのうち2ファイルを参照していた。

以下を `649b617` から復元し、Git blob の一致を確認した。

- `reviews/shift-distribution-check.json`
- `reviews/shift-distribution-check.log`
- `reviews/shift-response-changes.json`
- `reviews/shift-response-plan-2026-09-10.md`

これらは当時の記録であり、今回の版の検証結果ではない。現行文書にもその区別を
明記した。既存ソースを古い snapshot で上書きしてはいない。

## 検証

- Python 3.14: **4,139 passed, 82 skipped**、382.83秒。
- Python 3.10: **4,139 passed, 82 skipped**、402.50秒。
- 分岐を含む coverage: **85%**（必要基準80%）。
- Ruff lint/format、非 incremental mypy、依存整合性、Bandit、pip-audit、
  複雑度の hard limit: 全て通過。既知の脆弱性は検出されなかった。
- 独立した wheel/sdist のバイト再現性、内容、metadata 検査: 全て通過。
- インストール済み wheel: ARM/Intel 両 Python で checkout 外から import・
  version・validate・info・Newick/TSV 往復を確認。Intel 環境では SVG tip image
  の描画と、ARM でスキップした SVG 画像変換3テストも通過した。
- clean sdist: checkout 外で CLI・全コマンド handler・wiki example の
  **100テスト通過**。

| 環境 | Python | NumPy | pandas | SciPy |
|---|---|---|---|---|
| macOS ARM64・最低対応版 | 3.10.21 | 2.2.6 | 2.3.3 | 1.15.3 |
| macOS ARM64・最新対応版 | 3.14.7 | 2.5.3 | 3.0.5 | 1.18.1 |

両環境とも ETE4 4.4.0。`constraints-dev.txt` を用いて隔離環境を構築し、
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg`
を設定した。Python 3.14 では `python tools/check.py release`、
Python 3.10 では `python tools/check.py test` を実行した。

環境構築時、シェルに残っていた Intel 向け conda コンパイラ設定が ARM Python
の ETE4 拡張を壊すことを確認した。タスク専用環境で該当 compiler flags を外し、
ETE4 4.4.0 をソースから再構築して解消した。ユーザーの環境設定は変更していない。

IQ-TREE と R の kfl1ou は利用できず、それらを必要とする外部統合テストは
スキップ対象。ARM 環境の SVG rasterization 3テストも Intel 版 Cairo との
architecture 不一致によりスキップされ、Intel Python 環境で補完済み。
Windows/Linux の実行はローカルでは行わず、push 後の CI に委ねる。
認証情報の高確度パターンを送信対象の Python/JSON/YAML/TOML ファイルで検査し、
候補はなかった。最大の tracked blob は約13 MB で、100 MB 超のファイルはない。
