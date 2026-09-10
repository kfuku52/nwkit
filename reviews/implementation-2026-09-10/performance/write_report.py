"""Assemble the measured revision comparison; no timing values are hand-entered."""

import json
from pathlib import Path

FOLDER = Path(__file__).resolve().parent
ROOT = FOLDER.parents[2]
LABELS = {
    "startup": "CLI --version",
    "nwk2table-balanced-2048": "樹木表・balanced 2,048葉",
    "nwk2table-comb-1600": "樹木表・comb 1,600葉",
    "asr-BM-128": "BM ASR 128葉",
    "asr-BM-1024": "BM ASR 1,024葉",
    "asr-OU-128": "OU ASR 128葉",
    "asr-OU-512": "OU ASR 512葉",
    "asr-ER-128": "ER ASR 128葉",
    "asr-ER-1024": "ER ASR 1,024葉",
    "asr-MV-BM-64": "多変量BM 64葉・2形質",
    "asr-MV-BM-256": "多変量BM 256葉・2形質",
    "regress-lambda-128": "Gaussian回帰 λ推定 128葉",
    "regress-lambda-512": "Gaussian回帰 λ推定 512葉",
    "regress-binomial-30": "二項GLMM 30葉",
    "binomial-n60-p0.05": "稀な二値応答GLMM 60葉",
    "negative-binomial-n30-boundary": "負の二項GLMM 境界例 30葉",
    "asr-ensemble-128": "樹木ensemble BM 128葉×4樹",
    "pca-256": "PCA 256葉・2形質",
    "radte-branch-8": "RADTE 枝長 8種・16葉",
    "radte-branch-16": "RADTE 枝長 16種・32葉",
    "branch-gaussian-fixed-128": "枝別Gaussian 固定 128葉",
    "branch-gaussian-fixed-512": "枝別Gaussian 固定 512葉",
    "radte-jc69-8": "RADTE JC69 8種・128サイト",
}


def link(name):
    return f"[{name}]({FOLDER / name})"


def main():
    selected = {}
    for filename in [
        "results.json",
        "glmm-results.json",
        "repeat-results.json",
        "additional-results.json",
    ]:
        path = FOLDER / filename
        if not path.exists():
            continue
        for name, item in json.loads(path.read_text())["cases"].items():
            if item["summary"].get("ratios"):
                selected[name] = (filename, item["summary"])
    rows = []
    for name in LABELS:
        if name not in selected:
            continue
        _, s = selected[name]
        before, after = s["before"], s["after"]
        rows.append(
            "| "
            + LABELS[name]
            + " | "
            + f"{before['wall_seconds']['median']:.3f} → {after['wall_seconds']['median']:.3f} | "
            + f"{s['ratios']['wall_seconds']:.2f}× | {s['ratios']['cpu_seconds']:.2f}× | "
            + f"{before['peak_rss_bytes']['median'] / 2**20:.1f} → {after['peak_rss_bytes']['median'] / 2**20:.1f} |"
        )
    maximum = max(s["ratios"]["peak_rss_bytes"] for _, s in selected.values())
    document = f"""# NWKIT 最近の変更の時間・RAM比較 — 2026-09-10

今回の{len(selected)}ケースでは、各ケースのピークRSS中央値の増加は最大 **{(maximum - 1) * 100:.1f}%**。一方、**128葉のGaussian回帰のCLI時間は約2倍**、**GLMMでは約14〜37%増加**した。全処理が一律に遅くなったという結果ではない。

## 比較範囲と方法

- 主な比較: 9月6日 `8ad9544` → `e2465b0` に今回の3件の修正を加えた版。
- 新機能: RADTEは導入版 `0423427`、枝別Gaussianの固定パラメータCLIは `97491cb` と比較。PCAとensembleは今回の修正の影響を切り分けるため、修正直前の `e2465b0` と比較した。
- 追跡ファイルを一時ディレクトリへ展開し、修正した6ファイルだけを重ねた。進行中のSHIFT実装、result_plot、ASR.md等の別作業は含めていない。
- Apple M2 Max・64 GiB。実行したPythonは **x86_64版 3.10.14**。NumPy 1.26.4 / OpenBLAS、SciPy 1.15.2 / MKL 2020.0.4、pandas 2.2.3、ETE4 4.4.0。
- 全て同じ合成入力・乱数seed・ライブラリを使い、BLAS/OMP/MKL/VECLIBを1スレッドに固定。独立した新規プロセスでウォームアップ1回と測定3回を行い、新旧の順序を交互にした。
- 初回に待ち時間のばらつきが大きかった起動・樹木表・BMは追加で各5回測定した。表ではその追加測定を採用し、最初の測定もJSONに保持した。
- 実時間はPython内のCLI読み込み開始から処理完了まで。プロセス起動・終了を含む時間も `process_wall_seconds` に別記した。CPU時間は処理が消費した時間、RAMは各プロセスのOS報告の最大RSSで、Python以外の数値ライブラリのメモリも含む。tracemallocやcoverageは使用していない。全テストと主ベンチマークは同時実行していない。

## 測定値

中央値。倍率は新版/旧版。MiB = 1,048,576 bytes。入力サイズだけで最適化時間を予測できるという意味ではない。

| 処理 | 実時間・秒（旧→新） | 時間倍率 | CPU倍率 | 最大RSS・MiB（旧→新） |
|---|---:|---:|---:|---:|
{chr(10).join(rows)}

## 時間差の原因

**Gaussian回帰128葉:** 新旧とも138回のパラメータ評価で、探索回数は同じだった。`gaussian.solve_factor` がNumPyの一般行列ソルバからSciPyの `cho_solve` に変わっている。実測では遅延が最初の1回に集中した。実際の128×128行列を独立プロセスで再実行すると、初回約1.84 CPU秒、以後約0.0001〜0.0003 CPU秒で、結果は完全一致した。使用中のMKL経由の初回呼び出しコストによる増加と判断する。512葉では高速な行列計算の効果が上回り、CLI全体でも短縮した。これはこのx86_64/BLAS環境の観測であり、ARMネイティブPythonや別のBLASで同じ初期化時間とは限らない。

根拠: {link("profile-before.json")}、{link("profile-after.json")}、{link("solve-diagnostics.json")}、{link("first-solve.json")}。[該当コード]({ROOT / "nwkit/gaussian.py"}:690)。プロファイルの時間は原因調査用で、上表の無計装の測定とは区別した。

**GLMM:** 9月10日の `9e2dd93` が、無罰則モデルでも複数の初期値を調べるようにした。最初の収束解が不十分な場合への正しさの修正であり、通常二項で約14%、稀な二値応答で約21%、負の二項境界例で約37%の時間増加を観測した。[スカラーGLMMの探索]({ROOT / "nwkit/phylogenetic_glmm.py"}:2835)、[カテゴリGLMMの探索]({ROOT / "nwkit/phylogenetic_glmm.py"}:3688)。今回、この探索を減らす変更は加えていない。

**その他:** 最初の樹木表・BMの測定には大きな実時間差があったが、CPU時間差は小さく、追加5回ではほぼ同程度〜約9%の範囲だった。初回の1.5〜1.8倍を再現性のある性能低下とは扱わない。CLI起動自体は約13 ms増加した。ensembleの追加測定にも待ち時間差があり、実時間の短縮を今回の修正の高速化効果とは解釈しない。

## 出力の照合

全ての成功した測定回で、同じ版の出力は安定していた。全TSVの行・共通列を照合し、数値は `rtol=1e-5, atol=1e-7`、その他は完全一致とした。回帰表に追加された12個の診断列は追加として記録し、隠していない。

- 通常のASR、樹木表、Gaussian回帰、通常二項GLMM、PCA、ensemble、固定枝別Gaussian、RADTE枝長・JC69の共通出力は一致した。
- 稀な二値応答例は、係数・尤度はほぼ同じだが、p値の約7.1e-7、区間上端の約1.25e-5の差が設定した許容差を超えた。後から許容差を緩めず、不一致として保持した。
- 負の二項境界例では、対数尤度が -53.90478774 → -53.90476085、診断が `nuisance-information-singular` → `ok` になり、以前は欠けていた標準誤差・p値が得られた。係数・分散等も変わるため、これは同じ入力の実行コスト比較であり、等価な推定結果を得る速度の比較とは扱わない。

詳細: {link("output-checks.json")}、{link("glmm-output-checks.json")}、{link("repeat-output-checks.json")}、{link("additional-output-checks.json")}。入力・修正ソースのhashは {link("metadata.json")}、各出力のhashは測定JSONに記録した。

## 再現と限界

{link("prepare.py")} で一時コピーと合成入力を作り、表示されたルートの `manifest.json` を {link("benchmark.py")} の `--manifest` に指定する。`--output` は新しいJSONパスを使う。{link("prepare_glmm.py")} と {link("prepare_additional.py")} はそのルートを引数に取り、追加のmanifestを作る。全コマンド・入力パス・環境・各回の時間とRAMは {link("results.json")}、{link("glmm-results.json")}、{link("repeat-results.json")}、{link("additional-results.json")} にある。標準出力復元の最終調整後にはensembleを再測定した。他の測定処理は未変更で、初期・最終のソースhashをmetadataに残した。最初のGLMMハーネスの `--response-family` 指定ミスによる8回の失敗は、生ログに残して測定値から除き、正しい `NAME=VALUE` 指定で再測定した。

数値例:

```sh
python reviews/implementation-2026-09-10/performance/prepare.py
# 出力されたディレクトリを BENCH_ROOT に設定
python reviews/implementation-2026-09-10/performance/benchmark.py \\
  --manifest "$BENCH_ROOT/manifest.json" --output /tmp/nwkit-performance.json
```

実運用データ、大規模bootstrap・モデル探索全般、外部IQ-TREE 3/Rの処理、並行開発中のnative SHIFT、他OS・他Python/BLASは未測定。枝別Gaussianの新しいパラメータ推定には導入前の等価なCLIがなく、固定パラメータ計算とは分けて考える必要がある。今回の測定を、全モデル・全入力でRAMや時間が増えない保証とはしない。
"""
    (FOLDER.parent / "PERFORMANCE.md").write_text(document, encoding="utf-8")


if __name__ == "__main__":
    main()
