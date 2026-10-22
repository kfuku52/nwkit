# 第2回性能候補の本体採用

基準 `8664baa7e88afae7e384c0e4d8d16437846493a2` (0.43.56) に対し、
**7案を採用してそれぞれcommitしました。** [10案の独立評価](../performance-candidates-round2-2026-10-22/REPORT.md)
のうち6案に、固定rateへの制限を実装・検証したERを加えています。
WGD log product、threshold diagonal、Brownian factorの3案は保留のままです。

## 統合後の測定

下表は試作パッチを使わない、本体の操作全体の比較です。7回の交互比較を実施し、
各行のすべてのペアで改善、baseline/currentの生時間範囲も分離しました。
倍率を異なる操作にまたがって合成したり、全コマンドの効果とは解釈しません。

| 操作・条件 | 基準 ms | 本体 ms | 中央値倍率 | 出力 |
|---|---:|---:|---:|---|
| [λ fit、256 tip、profile CIあり](lambda.json) | 209.214 | 8.478 | 24.68× | 最大尺度付き差8.61e-9 |
| [相関PCA、512 tip・16 trait](pca.json) | 6.597 | 4.219 | 1.56× | 完全一致 |
| [Brownian REML PGLS、256 tip](pgls.json) | 22.107 | 6.611 | 3.34× | 完全一致 |
| [OU共分散作成、512 tip](covariance.json) | 96.455 | 12.720 | 7.58× | 完全一致 |
| [希少128家族、状態上限86](wgd.json) | 118.942 | 104.311 | 1.14× | 完全一致 |
| [Full-OU識別性診断、8 trait・32 tip・4096 moment](fullou.json) | 1781.160 | 246.007 | 7.24× | 完全一致 |
| [count map、32 tip・一意枝長・1 draw](maps.json) | 17.405 | 13.044 | 1.33× | 完全一致 |
| [history map、同条件・1 draw](histories.json) | 19.929 | 12.910 | 1.54× | 完全一致 |
| [ER全marginals、128 tip・64 state・固定rate](er.json) | 83.347 | 71.993 | 1.16× | 最大尺度付き差1.78e-15 |

Full-OUの倍率は有限差分・行列指数・SVDを含む識別性診断全体に対する値です。
自由MVOU fit全体が同じ倍率になるという主張ではありません。λ推定の既存の
探索範囲・11点grid・xatol・CI・検定は維持し、最適λの丸め差を許容する既存基準で
照合しました。固定profileは独立Gaussian密度とも整合します。
ERのlog likelihoodとrateは従来計算のまま、最終posteriorのみを加速しています。

[自由ER rate](er-free-control.json)、[自由λ PCA](pca-search-control.json)、
[100 map標本](maps-default-control.json)は、出力が完全一致しました。
速度倍率はそれぞれ0.98、1.02、1.02で生時間範囲が重なるため、速度差は主張しません。
[ゼロ係数の共分散](covariance-zero-control.json)も完全一致し、0.99倍で生範囲が
重なりました。安定化ループを維持し、試作にあった準備の再実行を除いています。

### メモリと限定条件

tracemallocによるピーク中央値で、PCAは8.67→6.30MB (約27%減)、
PGLSは4.22MBでほぼ不変、Full-OU診断は9.91MBでほぼ不変でした。
λ fitは2.70→3.17MB、OU共分散は2.52→3.08MB、希少128家族は1.69→2.73MBに
増えています。Gaussianの一時積は4MiB、WGDのbroadcast/logsumexp作業量は1MiBを
目安に分割し、大きな状態空間では従来の行処理を保ちます。これはRSSの上限や、
SciPy内部の全allocationを含む厳密なメモリ上限ではありません。

確率写像は事前の条件付きPoisson tail検証をcontext構築に引き渡し、一意のQ/枝長に
つき1回の準備にします。cacheは1要求内だけで破棄し、guard前の大配列確保を許しません。
新しい7回比較でcount・history・event timesが完全一致し、既存の並列/seed/rollback
テストも通過しました。default100 drawでは目立つ改善を主張しません。
ERは32状態以上・固定rate・厳密な対角定数＋非対角定数のPだけを対象にし、
1ULPだけ構造が異なる供給Pは汎用計算へ戻します。自由rateの探索には適用しません。

## 検証と再現

環境、warmup、測定順序、メモリ・同等性の定義は[PROTOCOL.md](PROTOCOL.md)、
固定archiveと初回統合本体のhashは[source-integrity.json](source-integrity.json)、
最終ソースhashは [VALIDATION.json](VALIDATION.json) に記録しました。
初回historyベンチはfixture結果の枝IDをノードとして索引するharnessの誤りで停止し、
枝IDによる索引に修正して当該操作を最初から7ペア測定しました。本体の問題ではありません。

各案の統合時に実行したcheck entrypointと結果は以下のログに残しました。

| 対象 | 成功テスト実行数 | ログ |
|---|---:|---|
| Full-OU、slow実optimizer含む | 10 | [fullou](checks/fullou.log) |
| PCA・数値不変性・CLI・画像 | 59 | [pca](checks/pca.log) |
| WGD count・独立参照・実bootstrap CLI | 106 | [wgd](checks/wgd.log) |
| Gaussian・evolution・数値不変性 | 96 | [covariance](checks/covariance.log) |
| signal・calibration・独立Gaussian・R参照 | 52 | [signal](checks/signal.log) |
| ASR・確率履歴・独立CTMC監査 | 85 | [stochastic](checks/stochastic.log) |
| ER制限・離散ASR・確率履歴 | 142 | [er](checks/er.log) |

重複を含む実行数で、uniqueテスト数の合計ではありません。各quickにはlint、format、
型検査を含みます。対象ファイル名の指定誤りでpytestが起動しなかった試行は成功に
含めず、実在する `tests/test_discrete_asr_models.py` を指定して改めて通過しました。
数値許容差・モデル・単位・SE・missing・tail・seed・work guardは変更していません。

追加の悪条件監査で、32 tip・2 block・対角gap 1e-8の例に尤度比差
2.164e-7を確認しました ([修正前12対照](conditioning-before.json))。初回releaseは
5254 pass・127 skipで通過しましたが、この境界を十分に検査していなかったため、
その結果を最終版の検証には使いません。

スペクトルを使う条件を `min_eigenvalue > n * sqrt(eps) * max_abs_eigenvalue`
に強め、悪条件では元のCholesky profileを維持しました。科学的rank、探索範囲や
許容差を変える変更ではありません。追加回帰は[5 pass](checks/conditioning-guard.log)、
[実装を使った36対照](conditioning-after.json)は最大差1.56e-9で基準内でした。
[再現スクリプト](conditioning_audit.py)は完全fit・CI・検定を比較します。
修正後のλベンチを最初から7ペア再測定し、24.68倍と確認しました。
[以前の25.25倍の生測定](lambda-before-conditioning-guard.json)も保存しています。

最終 `python tools/check.py release` は **5256成功・127 skip**、combined branch
coverage **85.96%**で通過しました。型・依存関係・Bandit・脆弱性・複雑性・
sdist/wheel内容・メタデータ・バイト再現性を含みます ([最終ログ](checks/release-final.log))。
生成wheelを別targetへ入れ、import先をassertした上で、sdist中のCLI/interfaceテストと
全26件の新規回帰を実行し、**154成功**を確認しました ([配布物consumer](package-consumer.json)、
[ログ](checks/package-consumer.log))。7対象moduleはcheckout・sdist・wheelでバイト一致しました。

検証SHA・元のremote SHA・最終ソースhash・artifact hash・実行環境・skip/unrunの
境界を [VALIDATION.json](VALIDATION.json) に記録しました。127件のskipは外部backend、
opt-in研究、platform条件等を含み、成功に読み替えていません。pushの最終SHAと
remote CIの結果はタスクの完了報告で示します。
母集団校正のopt-in研究や新規の外部R等による再実験は、この計算再利用変更の範囲では
実行していません。既存の独立参照とsmall real bootstrap回帰を維持しています。
