# NWKIT レビュー3件の修正と検証 — 2026-09-10

[前回レビュー](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/REVIEW.md)で再現した3件を修正した。[性能調査](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/PERFORMANCE.md)も完了した。

## 修正

| 問題 | 修正と検証 |
|---|---|
| ASRが `--regime-parameters` 入力を上書きする | 出力との衝突検査に当該入力を追加。主出力・model出力・tree出力について、同一パス・symlink・hardlinkの9ケースで入力保持を確認した。 |
| アンサンブル失敗後に新旧の結果が混在する | reference・ensemble・副出力を共通の出力トランザクションでstageする。tip集合の不一致とfit後の失敗で、既存ファイルを保持し、新規ファイルを残さないことを確認した。標準出力は計算完了まで一時ファイルに保留するため、出力サイズに比例するRAMバッファを追加しない。標準出力の書き込み・flush失敗もファイル側の復元対象とした。SVGと標準出力の成功経路も確認した。 |
| 200葉のcomb treeでPCAが停止する | 再帰的deepcopyを、保持対象ノードの反復的コピーに置き換えた。枝長・元のbranch IDを保持する。欠測なし／欠測を除外する2ケースでCLIが成功し、元の樹木から独立に計算した共分散・固有値と一致することを確認した。入力樹木のpropertyも変更しない。 |

実装: [asr.py](/Users/kf/repos/nwkit/nwkit/asr.py)、[asr_output.py](/Users/kf/repos/nwkit/nwkit/asr_output.py)、[asr_figure.py](/Users/kf/repos/nwkit/nwkit/asr_figure.py)、[pca.py](/Users/kf/repos/nwkit/nwkit/pca.py)。

回帰テスト: [test_asr_output.py](/Users/kf/repos/nwkit/tests/test_asr_output.py)、[test_pca.py](/Users/kf/repos/nwkit/tests/test_pca.py)。合計22ケースを追加した。Windowsでは権限が必要なsymlink作成ケースを既存テストの方針に合わせてスキップする。

既存の出力トランザクションは例外時の復元を提供するもので、プロセスクラッシュや複数ファイルの同時読み取りまで含む完全な原子性を保証するものではない。

## 検証結果

- 全体の `python tools/check.py full`: **3,772 passed / 79 skipped / 6 warnings**、終了コード0。branchを含むカバレッジ **85%**。
- その後の標準出力復元の最終調整とテスト整備について、ASR・ensemble・図・PCA・出力トランザクション関連を再実行: **136 passed**。Ruff lint/format、mypy（204ソースファイル）も再通過した。
- 最終ソースのBanditとmaintainability検査: 成功。Radon 3,151関数、平均6.75、最大50。既存baselineからの複雑度増加の警告はあるが、hard limitは全て通過した。
- 依存関係検査・脆弱性監査: 成功。
- 最終ソースの `python tools/check.py dist`: sdist/wheelの内容・再現性検査に成功。
- 性能調査: **23ケース**。最終stdout調整に影響されるensemble測定も再実行し、新旧の数値出力を照合した。

全体検査から最終調整までに変わったファイルと追加検証を明示しているため、カバレッジ85%を最終調整後の全テスト再測定値とは扱わない。最終ソースのhashは検証用コピーと作業ツリーで一致することを確認した。

検証は `e2465b0` の一時コピーに今回の変更だけを重ねて実行した。並行作業のSHIFT、README、result_plot等は保持し、混在させていない。配布物検査もこのコピーで実行し、作業ツリーの既存build/distを削除していない。コミット・pushは行っていない。

外部kfl1ou/IQ-TREE 3等を必要とする統合試験にはスキップがある。今回の実行環境はmacOS・Python 3.10.14で、他OSや他Pythonの検証ではない。

ログ: [全体検査](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/fixed-full-check.log)、[最終関連検査](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/fixed-final-targeted.log)、[配布物](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/fixed-dist-check.log)、[複雑度](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/fixed-maintainability.log)。

前回の `reproduce.py` / `results.json` は修正前の症状の記録として残している。修正後の確認には上記の回帰テストを使う。
