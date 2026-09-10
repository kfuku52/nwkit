# NWKIT 実装レビュー — 2026-09-10

> 追記: 下記3件は修正済み。[修正・検証記録](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/FIXES.md)と[性能調査](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/PERFORMANCE.md)を参照。以下は修正前のレビュー記録。

現行コードで、入力ファイルの上書き、失敗した実行の部分的な出力更新、PCAの有効な入力に対する停止を再現した。3件とも修正対象と判断する。ソースコードの修正は行っていない。

## 対象

- HEAD: `e2465b048f7529b7dad84cb6722278406bb49bbe`、version `0.43.15`。
- 主な差分範囲: `741836b^..e2465b0`（9月7日以降の38コミット）。`nwkit/` と `tests/` で235ファイル、44,859行追加・667行削除。変更された実装ファイルは146。
- 開始時から存在した `ASR.md` の未コミット変更と `reviews/scientific-2026-09-10/` の未追跡ファイルは保持した。以前の科学的レビューは旧HEADに対するものなので、その結論を現在の不具合として転載していない。
- 終盤に別作業による未追跡ファイル `nwkit/gaussian_whitening.py`、`nwkit/shift_native_fit.py`、`nwkit/shift_native_model.py` と対応テストの追加を検出した。これらは開始時のHEADとテスト収集に含まれず、本レビューの対象外。進行中の変更は保持した。
- 差分の整理、主要実装の静的レビュー、既存の全検査、新規の反例を組み合わせた。44,859行の全行に対する独立した数学的証明や、全プラットフォームでの実行を行ったという意味ではない。

## 指摘

### 1. [P1] `--regime-parameters` の入力ファイルが出力によって消える

対象: [nwkit/asr.py:2971](/Users/kf/repos/nwkit/nwkit/asr.py:2971)、`_validate_asr_output_paths`。

入力と出力の衝突検査の対象に `regime_parameters` が入っていない。そのため、BMSなどで `--regime-parameters parameters.tsv -o parameters.tsv` を指定すると、パラメータ表を読み込んだ後、警告・例外なしで同じファイルをASR結果表に置き換える。元の推定条件が失われる。他の入力ファイルはこの関数で保護しており、`asrcompare` は同じ `--regime-parameters` を保護している。

4葉・1レジーム・固定 `sigma2=1` で、元の `regime\tsigma2` ヘッダーが `branch_id\tparent\tnode_class...` に変わることをCLIから確認した。これは利用者が誤って出力先を入力に重ねた場合のデータ喪失であり、通常の出力先で勝手に入力を削除するという指摘ではない。

**修正案:** `regime_parameters` を保護対象に追加し、主出力・副出力、同一パス・symlink・hardlinkに対する拒否を確認する。この入力保護の欠落は現行で確認された残存不具合であり、今回の最新コミットだけで発生したと断定しない。

### 2. [P2] 樹木アンサンブルの失敗後に新旧の結果が混在する

対象: [nwkit/asr.py:3131](/Users/kf/repos/nwkit/nwkit/asr.py:3131)、[nwkit/asr_tree_ensemble.py:92](/Users/kf/repos/nwkit/nwkit/asr_tree_ensemble.py:92)。

通常のASR結果を直接書き出してから、アンサンブルの実際の樹木・tip集合・重み・各fitを検証する。BMなどの通常経路には出力一式のトランザクションがない。後半が失敗すると、主出力は今回の結果に更新済みなのに、既存の `--tree-ensemble-out` は前回のまま残る。

主樹を `A,B,C,D`、アンサンブルの樹を `A,B,C,E` にして、両出力に前回の内容を置いて実行した。`Every ensemble tree must contain exactly the reference tip set.` で失敗したにもかかわらず、主出力だけが更新された。解析失敗後に結果を参照すると、同じ実行に属さないファイルを組み合わせるおそれがある。

**修正案:** 入力集合などの事前検証を前倒しし、さらにreferenceとensembleを含む出力一式をstageして、全処理成功時にまとめて公開する。事前検証だけでは途中の最適化失敗・書き込み失敗を防げない。BRANCH-GAUSSIANと拡張stochastic mappingには既に同種の仕組みがある。

### 3. [P2] PCAが200葉の偏った樹形で `RecursionError` になる

対象: [nwkit/pca.py:67](/Users/kf/repos/nwkit/nwkit/pca.py:67)、`_retained_tree`。

保持する樹木を `copy(method="deepcopy")` で再帰的にコピーしている。この処理は欠測を落とす場合に限らず毎回行うため、欠測のない200葉のcomb treeでもPythonの標準再帰上限に到達して停止する。

CLIから再現し、同じ樹木・同じ2形質をコピー処理を経由せずPCAの数値関数に渡すと `status=ok`、固有値 `[422.0261174239691, 0.9968698249142601]` を得た。データのrank不足や共分散の非正定値性ではない。限界の正確な葉数は樹形とPython環境に依存する。

**修正案:** 枝長と元のbranch IDを保持する反復的なコピーを使う。DTTの `_prune_to_crown` は既に同じ問題を避けている。再帰上限を引き上げるだけの対処は避ける。

## 確認した範囲

| 領域 | 重点確認 |
|---|---|
| BRANCH-GAUSSIAN | 固定・推定パラメータの受け渡し、root条件、共有グループ・境界、中心化、識別性検査、最適化の確認、出力トランザクション |
| ASR拡張 | モデル平均とtree ensemble、観測誤差・replicateの尤度定数、CVでのholdout除去、再fitの固定/自由パラメータ契約、posterior/予測/bootstrap出力 |
| 多変量・個体ASR | covarianceの正規化・復元、full OUの固定行列、個体と種の推定対象、欠測とSEの扱い |
| 回帰 | event-averageと共通係数の分離、censoring機構の再生成、正則化推論の契約、nested CVとtraining内標準化、条件付き予測 |
| RADTE | rooted/unrooted枝長の対応、codon/DNA/AA likelihood、gamma categories、quadratic近似検証、input ensembleと区間の扱い |
| SHIFT・THRESHOLD | backend能力検証、推論上の制限の表示、rank/split診断と監視対象、既存テスト |
| PCA・signal・DTT・従来コマンド | 単位の扱い、欠測時の樹木処理、出力保護、集合/樹木操作の差分、配布パッケージ |

この表はレビューの焦点を示し、各領域に未知の問題がないことを保証しない。上の3件以外について、この作業で報告できる再現済みの新しい不具合は得ていない。

## 実行記録

環境: macOS、Python 3.10.14、NumPy 1.26.4、SciPy 1.15.2、pandas 2.2.3、ETE4 4.4.0。

実行コマンド:

```sh
python tools/check.py full
python tools/check.py dist
python -m ruff check reviews/implementation-2026-09-10/reproduce.py
python reviews/implementation-2026-09-10/reproduce.py
```

`dist` はHEADの追跡ファイルを別ディレクトリに展開して実行した。元の `build/`、`dist/`、`direct-dist/` を削除せず、元のファイル権限と `SOURCE_DATE_EPOCH` を保持した。最初の展開で生じた権限差によるwheelのバイト差は、元の権限に揃えて再実行し解消した。パッケージの内容差はなかった。

検査結果:

- `full`: 終了コード0。**3,752 passed、79 skipped、6 warnings**、テスト実行27分56秒。
- Ruff lint/format、mypy（203ソースファイル）、`pip check`、Bandit、`pip-audit` は成功。既知の脆弱性は検出されなかった。
- branch coverageを有効にした総合カバレッジは、レビュー対象で **85%**。元のログは実行中に追加された未追跡3モジュールを0%として含み **84%** だったため、同じ測定データから当該3ファイルだけを除いて再集計した。テストの再実行はしていない。
- maintainabilityのhard limitは全て通過。Radon集計は3,194関数・平均6.76・最大50。既存baselineからの複雑度増加の警告が残る。この終盤の集計には並行作業のファイルが含まれ得る。
- `dist`: sdist/wheelの内容・再現性検査に成功。
- 追加の3件の再現assertと、再現スクリプトのRuff検査は成功。

79件のスキップがあり、全外部連携を実行したわけではない。特にRの `kfl1ou` は未インストール、`NWKIT_TEST_RSCRIPT` は未設定、IQ-TREE 3は利用できず、これらを必要とする実機統合は未検証。6件のwarningはrequestsの依存バージョン警告、テスト内pandas concatの将来変更、数値最適化時の差分・境界警告だった。監査時にはキャッシュ読み取り警告もあったが、監査は完了した。

完全なログ: [full-check.log](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/full-check.log)、[dist-check.log](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/dist-check.log)。

新しい大規模な型I誤差・被覆率実験、代表的条件での新旧性能比較、Windows/Linux、Python 3.11〜3.14での実行、図の全ページの目視検査は行っていない。テストの成功を科学的な校正や全入力での正しさの証明とは扱わない。

## 再現資料

- [reproduce.py](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/reproduce.py): 全て一時ディレクトリ内のダミーファイルで再現する。リポジトリの入力データを上書きしない。現在の3不具合が再現することをassertする診断スクリプトであり、修正後に通る回帰テストではない。
- [results.json](/Users/kf/repos/nwkit/reviews/implementation-2026-09-10/results.json): 実行結果、HEAD、対象ソースのSHA-256。

優先順は、入力データの保護、出力一式の整合性、PCAの樹形制約の解消。
