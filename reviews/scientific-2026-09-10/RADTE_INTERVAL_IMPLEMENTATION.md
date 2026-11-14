# RADTE区間推定の対応結果

初期計画に沿って実装し、検証結果により採用範囲を限定した。既定の
`--uncertainty none` とstudentizedの数式は維持した。commit/pushはしていない。
参照元 `/Users/kf/repos/nwkit` は編集せず、割当worktreeだけを変更した。

## 実装したもの

1. **marginalのゼロ分散境界**：最適化変数をlog SDから非負分散に変更した。
   SDの数値下限を小さくする代わりに、分散0をモデルに含め、root pairの
   非線形性を含む右微分を導出・実装した。科学的なゼロ分散境界と数値上の
   上限・平均の限界を区別する。小さい正分散で有限差分が負分散側へ出ない
   ようにし、変数の単位変更だけで「境界」と判定しないようにした。
2. **限定された厳密区間**：通常CLIに `--uncertainty exact-log-duration`
   を追加した。枝長入力、自由年齢1個、変動するdurationが共通の
   `age-offset` の正の倍数となる場合だけ利用できる。未知分散では厳密な
   GLSのt対比、既知の正SDでは正規対比を使う。適用不能な樹や配列入力は
   拒否する。点推定値は変えず、既存の区間列・status・manifestで出力する。
3. **一般校正profileの研究用実装**：制約下のnuisance再推定、rates→配列
   の二段階生成、marginal尤度比、失敗を捨てないMonte Carlo上限、全gridを
   見る区間候補を実装した。ただし下記の失敗率により、通常CLIの選択肢には
   採用しなかった。単にMAPへfallbackしたり近似検査を外したりしていない。
4. **検証runner**：樹形、生成/推定rho、生成/推定置換モデル、calibration幅、
   推定方式、枝長のみの実験、複数target、optimizer seed、開発/検証の区別を
   追加した。ある区間方式だけの例外で他方式の結果を消さず、返却率R/N、
   条件付き被覆C/R、正しい区間の返却割合C/N、理由別失敗、二項区間を保存する。
5. **独立参照と証拠**：主実装と別のdense GLS参照、dense共分散/expm生成器、
   family-level CSV、protocol、source hash、開発pilotを保存した。

## 検証結果と採否

- 枝長入力の3条件・各1,000独立試行で、厳密経路は全例に区間を返し、
  真値の包含は960、946、942例。独立参照との最大端点差は正規化時間で
  `3.11e-15`。Laplaceは同じ試行で809、803、809例だった。
- ただし事前に定めた同時片側下限≥.93の基準は、厳密経路でも後2条件が
  .9287、.9242で未達。**実験の合格とはしていない**。追加seedや乗数変更で
  通過させず、限定モデルの数理的な厳密性・数値一致と、実験上の採用基準を
  分けて記録した。
- 配列の低SD・marginal限定20例では、17例が二次近似のexactチェックで停止。
  残る3例は分散0であり、曲率区間は利用不能だった。これは元のauto解析の
  成績ではなく、marginal限定の開発pilotである。
- その1例に対する真の年齢での二段階帰無生成19回でも12回が近似チェックに
  失敗した。したがって一般的な校正profileの利用可能性を満たさず、研究用に
  留めた。この結果を校正済みP値や被覆率と呼んでいない。

詳細と再現コマンドは
[保存証拠](../../examples/radte/interval-boundary-validation/README.md) にある。
記録後のstatus分類やstrict-clockガードの追加も明記し、保存版のhashを
現在版のhashと偽っていない。

## 残る科学的課題

一般の配列データ、複数/内部重複、広いcalibration、置換・clock誤指定で
過小被覆と利用不能率が解消したとは言えない。現段階で一般向け95%方式を
推奨する判断は見送る。次の根本課題は、現在停止する領域でも誤差を管理できる
sequence marginal計算と、nuisance plug-in/grid近似の検証である。
独立検証の失敗を根拠なく許容する変更は行わない。

topology誤り、gamma/codon誤指定、複数重複生成、全計画セルの大規模な独立検証は
未実施。MCMCTreeの出力差を正確性とは扱わず、今回も正解基準として使っていない。

## ソフトウェア検証

`tools/check.py quick` でformat・Ruff・mypyを通し、RADTE本体、marginal、
profile、studentized、配列/codon、species/ensemble、IQ-TREE/PAML、CLI契約の
対象テストは **234 passed / 27 skipped / 0 failed**。追加のstrict-clock
区間拒否テストも **1 passed**。利用可能なPAML実行の小規模契約テストは含むが、
MCMCTreeによる被覆研究ではない。依存関係のwarningとSLSQPの一時的な境界clip
warningが残った。maintainabilityのhard limitと`git diff --check`も通過した。
全repositoryのfull/配布チェックは実行していない。
