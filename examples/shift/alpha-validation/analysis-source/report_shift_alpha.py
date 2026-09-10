"""Create a canonical portable-report input from audited OU simulation results."""

import argparse
import gzip
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path


def build(folder):
    summaries = json.loads((folder / "summary.json").read_text())
    paired = json.loads((folder / "paired.json").read_text())
    audit = json.loads((folder / "audit.json").read_text())
    with gzip.open(folder / "records.jsonl.gz", "rt") as handle:
        records = [json.loads(line) for line in handle]
    refit_gaps = [
        abs(r["same_model_refit_score_gap"])
        for r in records
        if r["status"] == "completed"
        and r["method"] == "two_stage"
        and r["same_model_refit_score_gap"] is not None
    ]
    diagnostic_path = folder / "null-profile-diagnostic.json"
    diagnostic = (
        json.loads(diagnostic_path.read_text()) if diagnostic_path.exists() else None
    )
    if audit["point_models_audit_failed"]:
        raise ValueError(
            "Resolve selected-model audit failures before publishing this report"
        )
    blocks, charts, tables, datasets = [], [], [], {}
    title = "OU シフト選択：独立反復と α 下限への感度"

    def markdown(key, body, sourced=True):
        blocks.append(
            {
                "id": key,
                "type": "markdown",
                "body": body,
                **({"sourceId": "simulation"} if sourced else {}),
            }
        )

    def table(key, title, rows, columns):
        datasets[key] = rows
        tables.append(
            {
                "id": key,
                "title": title,
                "dataset": key,
                "sourceId": "simulation",
                "defaultSort": {"field": columns[0][0], "direction": "asc"},
                "columns": [
                    {"field": field, "label": label, "type": "text"}
                    for field, label in columns
                ],
            }
        )
        blocks.append({"id": key + "-block", "type": "table", "tableId": key})

    def chart(key, title, rows, subtitle):
        datasets[key] = [
            {**row, "floor": label, "rate": row[floor]}
            for row in rows
            for floor, label in (("small", "αH 下限 10⁻⁷"), ("raised", "αH 下限 0.1"))
        ]
        charts.append(
            {
                "id": key,
                "title": title,
                "subtitle": subtitle,
                "showDescription": True,
                "type": "bar",
                "dataset": key,
                "sourceId": "simulation",
                "layout": "full",
                "encodings": {
                    "x": {
                        "field": "label",
                        "type": "nominal",
                        "label": "選択基準 / 探索",
                    },
                    "y": {
                        "field": "rate",
                        "type": "quantitative",
                        "label": "割合",
                        "format": "percent",
                    },
                    "color": {"field": "floor", "type": "nominal", "label": "α 下限"},
                },
                "palette": {"kind": "categorical", "name": "blue-yellow"},
                "valueFormat": "percent",
                "settings": {"groupMode": "grouped", "sort": "none"},
                "labels": {"values": "all"},
                "referenceLines": [
                    {"axis": "y", "value": 1, "label": "100%", "color": "neutral"}
                ],
            }
        )
        blocks.append({"id": key + "-block", "type": "chart", "chartId": key})

    def selected(**conditions):
        return [r for r in summaries if all(r[k] == v for k, v in conditions.items())]

    def rate_text(row, metric):
        entry = row[metric]
        if entry["rate"] is None:
            return "推定不能"
        lo, hi = entry["wilson_95"]
        return f"{entry['count']}/{entry['denominator']} ({entry['rate']:.0%}; 95% CI {lo:.0%}–{hi:.0%})"

    markdown("title", "# " + title, False)
    null = selected(family="primary", scenario="null", method="joint")
    null_ranges = {
        criterion: [r["any_shift"]["rate"] for r in null if r["criterion"] == criterion]
        for criterion in ("BIC", "pBIC")
    }
    markdown(
        "summary",
        "## 技術要約\n\n"
        f"**独立な {audit['datasets_planned']} データセットで比較を完了しました。** "
        "主解析は 8 末端、4 種類の真値、2 種類の根モデルで各 50 反復です。各データに BIC/pBIC、二段階/joint、2 通りの α 下限を適用しました。\n\n"
        + "主解析の無シフト真値に対する joint の偽陽性率は、根と下限によって "
        + "; ".join(f"{c}: {min(v):.0%}–{max(v):.0%}" for c, v in null_ranges.items())
        + " でした。これはこの条件での選択頻度であり、名目有意水準を保証する検定ではありません。"
        + f"選択モデル {audit['point_models_completed']:,} 件のうち独立計算による監査不一致は {audit['point_models_audit_failed']} 件でした。"
        + "正しい尤度の計算と、適切なモデル選択は別の検証課題です。",
    )
    markdown(
        "definitions",
        "## 偽陽性・回復率の定義\n\n"
        "偽陽性は、真のシフトがないデータで、親と異なる最適値を持つ実効シフトを 1 個以上選ぶことです。"
        "背景と統合されて消えた枝は数えません。回復率は、末端を同じ最適値でまとめた分割が真値と完全一致する割合です。"
        "これは枝位置の完全一致と異なります。RMSE は観測値への当てはまりでなく、生成過程の真の末端平均に対する誤差です。\n\n"
        "割合の分母は各条件で完了した反復数。失敗数と全試行を分母にした保守的な上下限は集計データに残しています。"
        "各棒は同じデータから得た対応のある結果で、棒どうしを独立標本とは扱いません。95% 区間は Wilson 区間です。",
    )
    for root, name in (("OUfixedRoot", "固定根"), ("OUrandomRoot", "定常ランダム根")):
        rows = selected(family="primary", scenario="null", root_model=root)
        evidence = []
        for criterion in ("BIC", "pBIC"):
            for method in ("two_stage", "joint"):
                values = [
                    r
                    for r in rows
                    if r["criterion"] == criterion and r["method"] == method
                ]
                entry = {
                    "label": f"{criterion} / {'二段階' if method == 'two_stage' else 'joint'}",
                    "root_model": root,
                }
                for value in values:
                    entry[value["floor_id"]] = value["any_shift"]["rate"]
                    entry[value["floor_id"] + "_details"] = rate_text(
                        value, "any_shift"
                    )
                    entry[value["floor_id"] + "_count"] = value["any_shift"]["count"]
                    entry[value["floor_id"] + "_n"] = value["completed"]
                evidence.append(entry)
        bullets = [
            f"{r['criterion']}・下限 {r['floor_id']}: {rate_text(r, 'any_shift')}"
            for r in rows
            if r["method"] == "joint"
        ]
        markdown(
            "null-" + root,
            f"## {name}：無シフト真値での選択\n\n"
            "棒の高さが高いほど、存在しない最適値の変化を選びやすいことを示します。"
            "下限の引き上げは推定範囲そのものを変えるため、結果を見て採用する調整ではなく感度分析として読みます。\n\n"
            "joint の結果：\n\n" + "\n".join("- " + x for x in bullets),
        )
        chart(
            "false-positive-" + root,
            name + "：偽陽性率",
            evidence,
            "8 末端・無シフト真値・各条件 50 反復。区間と分母は直前の本文とデータに保存。",
        )
    for root, name in (("OUfixedRoot", "固定根"), ("OUrandomRoot", "定常ランダム根")):
        rows = selected(family="primary", scenario="convergent", root_model=root)
        evidence = []
        for criterion in ("BIC", "pBIC"):
            for method in ("two_stage", "joint"):
                entry = {
                    "label": f"{criterion} / {'二段階' if method == 'two_stage' else 'joint'}",
                    "root_model": root,
                }
                for r in rows:
                    if r["criterion"] == criterion and r["method"] == method:
                        entry[r["floor_id"]] = r["shared_recovered"]["rate"]
                        entry[r["floor_id"] + "_details"] = rate_text(
                            r, "shared_recovered"
                        )
                        entry[r["floor_id"] + "_rmse"] = r["mean_rmse"]
                evidence.append(entry)
        markdown(
            "recovery-" + root,
            f"## {name}：真の収斂を回復できるか\n\n"
            "この条件では、離れた 2 クレードが同じ最適値を持ちます。棒が高いほど共有最適値の末端分割を正しく回復しています。"
            "探索を広げても、選択基準が真の分割を優先するとは限りません。\n\n"
            + "\n".join(
                f"- joint {r['criterion']}・下限 {r['floor_id']}: {rate_text(r, 'shared_recovered')}"
                for r in rows
                if r["method"] == "joint"
            )
            + "\n\n境界到達（joint、下限 10⁻⁷）："
            + " / ".join(
                f"{r['criterion']} は下限 {r['lower_boundary']['count']}/{r['completed']}、上限 {r['upper_boundary']['count']}/{r['completed']}"
                for r in rows
                if r["method"] == "joint" and r["floor_id"] == "small"
            )
            + "。上限への集中もあるため、この下限比較をパラメータ範囲全体に対する頑健性の証明とは扱いません。",
        )
        chart(
            "convergence-" + root,
            name + "：共有最適値の回復率",
            evidence,
            "8 末端・収斂真値・各条件 50 反復。完全な末端分割一致を評価。",
        )
    markdown(
        "paired",
        "## 下限による変化を同じデータの中で比較\n\n"
        "次の表は αH の下限を 10⁻⁷ から 0.1 に変えたときの分割変更数です。"
        "RMSE 差は引き上げ後−引き上げ前で、負なら真の平均に近づきます。"
        "分割が変わっても必ずしも正解に近づくとは限りません。主解析の各行は 50 組です。",
    )
    changes = [
        r for r in paired if r["family"] == "primary" and r["dimension"] == "floor_id"
    ]
    table(
        "floor-pairs",
        "下限変更の対応比較",
        [
            {
                "condition": f"{r['scenario']} / {r['root_model']}",
                "method": r["method"],
                "criterion": r["criterion"],
                "changed": f"{r['partition_changed']['count']}/{r['completed']}",
                "recovery_gain": r["shared_recovered"]["right_only"],
                "recovery_loss": r["shared_recovered"]["left_only"],
                "rmse_delta": f"{r['mean_rmse_right_minus_left']:.5f}",
            }
            for r in changes
        ],
        [
            ("condition", "真値 / 根"),
            ("method", "探索"),
            ("criterion", "基準"),
            ("changed", "分割変更"),
            ("recovery_gain", "回復へ改善"),
            ("recovery_loss", "回復を喪失"),
            ("rmse_delta", "RMSE 差"),
        ],
    )
    markdown(
        "method-pairs",
        "## 探索を広げた効果と基準を変えた効果\n\n"
        "次の表は各真値・根・下限での対応比較です。『改善』は左では誤り、右では正しい分割になった件数、"
        "『喪失』はその逆です。二段階→joint と BIC→pBIC を別の比較として示します。"
        "これにより、候補の最小スコアを下げることと真の分割を回復することを区別できます。"
        "RMSE 差は右−左です。両基準間でスコアそのものは比較しません。",
    )
    method_rows = [
        r
        for r in paired
        if r["family"] == "primary" and r["dimension"] in ("method", "criterion")
    ]
    table(
        "method-criterion-pairs",
        "探索法・選択基準の対応比較",
        [
            {
                "condition": f"{r['scenario']} / {r['root_model']} / {r['floor_id']}",
                "contrast": f"{r['left']} → {r['right']}",
                "fixed": r.get("criterion", r.get("method")),
                "changed": f"{r['partition_changed']['count']}/{r['completed']}",
                "gain": r["shared_recovered"]["right_only"],
                "loss": r["shared_recovered"]["left_only"],
                "rmse": f"{r['mean_rmse_right_minus_left']:.5f}",
                "score_lower": (
                    f"{r['search_score_comparison']['joint_lower']}/{r['search_score_comparison']['eligible']}"
                    if "search_score_comparison" in r
                    else "比較対象外"
                ),
            }
            for r in method_rows
        ],
        [
            ("condition", "真値 / 根 / 下限"),
            ("contrast", "対応比較"),
            ("fixed", "共通設定"),
            ("changed", "分割変更"),
            ("gain", "回復へ改善"),
            ("loss", "回復を喪失"),
            ("rmse", "RMSE 差"),
            ("score_lower", "joint スコア低下 / 比較可能組"),
        ],
    )
    extension_findings = []
    for family, label in (
        ("weak_pull", "弱い復元力"),
        ("known_error", "既知の観測誤差"),
        ("sixteen_tips", "16 末端"),
    ):
        clauses = []
        for criterion in ("BIC", "pBIC"):
            rates = [
                r["any_shift"]["rate"]
                for r in selected(
                    family=family,
                    scenario="null",
                    method="joint",
                    floor_id="small",
                    criterion=criterion,
                )
            ]
            if rates:
                clauses.append(
                    f"{criterion} の偽陽性率 {min(rates):.0%}–{max(rates):.0%}"
                )
        if clauses:
            extension_findings.append(f"- {label}：" + "、".join(clauses) + "。")
    markdown(
        "all-results",
        "## 弱い復元力・観測誤差・16 末端での補助検証\n\n"
        "補助条件は各 10 反復の探索的検証です。主解析と混ぜて単一の率に集約せず、以下の全条件表で確認できます。"
        "10 反復では区間が広く、根・樹形・標本数一般への外挿はできません。境界頻度はパラメータの弱い識別性を読むために併記します。\n\n"
        "joint・下限 10⁻⁷ の無シフト真値での範囲（2 種類の根、各 10 反復）：\n\n"
        + "\n".join(extension_findings)
        + "\n\n補助条件でも過剰選択は残り、16 末端に増やすことだけで一貫して解消したとはいえません。",
    )
    table(
        "all-cells",
        "条件別の回復率・誤差・境界頻度",
        [
            {
                "condition": f"{r['family']} / {r['scenario']} / {r['root_model']}",
                "fit": f"{r['criterion']} / {r['method']} / {r['floor_id']}",
                "recovery": rate_text(r, "shared_recovered"),
                "any_shift": rate_text(r, "any_shift"),
                "rmse": f"{r['mean_rmse']:.5f}" if r["mean_rmse"] is not None else "NA",
                "lower": f"{r['lower_boundary']['count']}/{r['completed']}",
                "upper": f"{r['upper_boundary']['count']}/{r['completed']}",
                "failed": r["failed"],
            }
            for r in summaries
        ],
        [
            ("condition", "条件"),
            ("fit", "選択方法 / 下限"),
            ("recovery", "分割回復 (95% CI)"),
            ("any_shift", "実効シフトあり (95% CI)"),
            ("rmse", "平均 RMSE"),
            ("lower", "下限到達"),
            ("upper", "上限到達"),
            ("failed", "失敗"),
        ],
    )
    markdown(
        "design",
        "## 生成モデルと比較範囲\n\n"
        "独立な枝ごとの OU 革新でデータを生成しました。主解析の αH=2.1、σ²H=0.25、非ゼロの最適値は ±2。"
        "補助条件は αH=0.2、既知の観測標準誤差 0.2、または 16 末端です。H は樹高です。"
        "真値は無シフト・単一・異なる 2 最適値・共有する 2 最適値の 4 種。主解析は各真値×根で 50 反復、補助条件は無シフトと収斂のみ各 10 反復です。\n\n"
        "探索の上限は 2 シフト。二段階は非制約モデルの枝配置を全探索してから後退統合、joint は識別可能な枝配置と最適値の等値制約を同時列挙します。"
        "候補数は 8 末端で 195、16 末端で 899。両者とも αH 上限 10、開始値 1、下限 10⁻⁷ または 0.1。"
        "以前の既定上限は生成値を含まなかったため、この検証では明示的な共通上限を使いました。以前の実験との単純な差分比較はしません。\n\n"
        "seed は 20260911 を基に条件と反復から事前生成し、全入力を推定前に固定しました。過去の 12 データの再利用ではありません。"
        "補正済み・未リリースの kfl1ou 3.0.9 を隔離ライブラリで使用。通常の同バージョン番号のリリース版とは異なります。",
    )
    markdown(
        "validation",
        "## 独立監査と残る数値的限界\n\n"
        f"選択モデルは {audit['point_models_completed']:,}/{audit['point_models_attempted']:,} 件完了。"
        f"joint の候補試行は {audit['candidate_attempted']:,} 件、候補失敗は {audit['candidate_failed']} 件。"
        f"Python で再計算した平均の最大絶対差は {audit['max_mean_error']:.3g}、尤度は {audit['max_likelihood_error']:.3g}、"
        f"BIC は {audit['max_bic_error']:.3g} でした。監査許容差は 10⁻⁶ です。\n\n"
        "候補を全列挙しても連続パラメータの大域最適性は保証されません。境界での pBIC は α の小さな数値差にも敏感で、"
        "同一モデルの再最適化差を探索の改善と混同しないよう元記録に分離して保存しました。"
        f"同一モデルのスコア差は最大 {max(refit_gaps):.6g} で、差が 10⁻⁴ を超えた {sum(gap > 1e-4 for gap in refit_gaps)} 組を探索スコアの比較から除外しました。"
        "候補失敗があるデータについては、全候補での最小値を証明したとは扱いません。\n\n"
        "Wilson 区間は各率の標本誤差を表します。条件が多いため同時被覆率ではありません。"
        "この実験は平衡二分樹、単一形質、最大 2 シフトに限られ、ブートストラップ被覆率やモデル平均の校正は評価していません。",
    )
    if diagnostic:
        markdown(
            "null-optimization",
            "## 無シフトモデルの単純な最適化失敗では説明できない\n\n"
            f"高い偽陽性率を確認した後の追加診断として、主解析の無シフト真値を {len(diagnostic['records'])} 条件で独立に調べました。"
            "切片と分散を解析的にプロファイルし、α の対数間隔 161 点と各格子局所最大の周辺を精査しました。"
            f"保存された無シフトモデルからの尤度改善の最大値は {diagnostic['maximum_log_likelihood_improvement']:.3g} で、丸め誤差程度でした。"
            "この点検は元の選択結果を変更せず、事前計画した主解析とも区別しています。"
            "連続大域最適性の厳密な証明ではなく、シフトを含む全候補の追加多点最適化も実施していません。",
        )
    markdown(
        "next",
        "## 次に検証すべきこと\n\n"
        "- α 下限と境界付近の pBIC の振る舞いを先に扱い、今回の結果だけでモデル平均の重みを確率として採用しない。\n"
        "- 有限標本での過剰選択に対して、独立な検証データを用いた基準の校正を検討する。今回の検証データで調整したら、新しい seed で再検証する。\n"
        "- 実データの樹形・末端数・誤差条件に合わせた模擬実験で再現性を確認する。",
    )
    markdown(
        "questions",
        "## 残る問い\n\n"
        "非平衡樹や多形質で同じ傾向が残るか。小さい α での最適値の非識別性に、選択基準をどう対応させるべきか。"
        "真の平均の予測を目的とする場合、枝・分割の完全回復より良い評価法があるか。これらは今回の結果だけでは決まりません。",
    )
    # The portable renderer requires actual SQL provenance. Materialize the
    # reviewed Python aggregates, execute this read, and use its returned rows.
    (folder / "report-data.json").write_text(
        json.dumps(datasets, ensure_ascii=False) + "\n"
    )
    query = "SELECT dataset, position, payload FROM reviewed_rows ORDER BY dataset, position"
    with sqlite3.connect(":memory:") as connection:
        connection.execute(
            "CREATE TABLE reviewed_rows (dataset TEXT, position INTEGER, payload TEXT)"
        )
        connection.executemany(
            "INSERT INTO reviewed_rows VALUES (?, ?, ?)",
            [
                (key, i, json.dumps(row, ensure_ascii=False))
                for key, rows in datasets.items()
                for i, row in enumerate(rows)
            ],
        )
        datasets = {}
        for key, _, payload in connection.execute(query):
            datasets.setdefault(key, []).append(json.loads(payload))
    source = {
        "id": "simulation",
        "label": "固定プロトコルによる独立 OU 模擬実験",
        "path": "report-data.json",
        "query": {
            "engine": "SQLite (in-memory reviewed aggregates from Python and R)",
            "language": "sql",
            "sql": query,
            "description": "Paired simulation; aggregate completed fits within each frozen cell; preserve all failures.",
            "tables_used": [
                "reviewed_rows",
                "summary.json",
                "paired.json",
                "audit.json",
                "records.jsonl.gz",
                "candidate-ledger.jsonl.gz",
                "protocol.json",
            ],
            "filters": [
                "520 independent datasets; same datasets across all eight fit settings; no outcome-based stopping"
            ],
            "metric_definitions": {
                "any_shift": "at least one effective optimum change",
                "shared_recovered": "exact label-invariant tip partition",
                "rmse": "root mean squared error against generating tip means",
            },
        },
    }
    stamp = datetime.now(timezone.utc).isoformat()
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": title,
            "generatedAt": stamp,
            "blocks": blocks,
            "charts": charts,
            "tables": tables,
            "cards": [],
            "sources": [source],
        },
        "snapshot": {
            "version": 1,
            "generatedAt": stamp,
            "status": "ready",
            "datasets": datasets,
            "accessIssues": [],
        },
        "sources": [source],
    }
    (folder / "artifact.json").write_text(
        json.dumps(artifact, indent=2, ensure_ascii=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    build(parser.parse_args().folder)
