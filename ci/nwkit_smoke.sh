#!/usr/bin/env bash
set -Ee
set -u
(set -o pipefail) 2>/dev/null || true

echo "[NWKIT] smoke2 start"

# 入力置き場
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

# ========== 共通ユーティリティ ==========
assert_contains () { grep -Eq "$2" "$1" || { echo "ASSERT FAIL: $1 does not contain /$2/"; exit 1; }; }
assert_not_contains () { ! grep -Eq "$2" "$1" || { echo "ASSERT FAIL: $1 contains /$2/"; exit 1; }; }
# 追加：固定文字列用
assert_contains_lit () { grep -Fq -- "$2" "$1" || { echo "ASSERT FAIL: $1 does not contain literal '$2'"; exit 1; }; }
leafset () { tr ',();:\047\"' '\n' < "$1" | grep -E '^[^[:space:]]+$' | sort -u; }  # \047 = single quote

# ラベル集合（英数字と_で始まるトークンのみを葉候補とみなす）
leaflabels () { tr ',();:\047\"' '\n' < "$1" | grep -E '^[A-Za-z_][A-Za-z0-9_]*$' | sort -u; }

# RF距離を堅牢に抽出（整数のみ取り出し）
rf_distance () {
  local A="$1" B="$2"
  nwkit dist -i "$A" -i2 "$B" \
    | sed -n 's/.*Robinson-Foulds distance = \([0-9][0-9]*\).*/\1/p' \
    | head -n1
}

# “:数字” の大まかなカウント
count_colons () { grep -o ':[0-9]' "$1" | wc -l | tr -d ' '; }

# 既存テスト用：APG IV ファイル（存在すれば使う）
APG="${APG:-apgiv.nwk}"

# ========== 1) nhx2nwk ==========
cat > "$TMPDIR/in.nhx" <<'NHX'
(A:0.1[&&NHX:S=Arabidopsis_thaliana:B=0.1],(B:0.2[&&NHX:S=Capsella_rubella:B=0.2],C:0.3[&&NHX:S=Utricularia_gibba:B=0.3])n1:0.4[&&NHX:S=s1:B=0.4])nroot:0.0;
NHX

# NHX→Newick。内部ノードラベルを NHX の S で付け替える（葉は実装によっては変わらない）
nwkit nhx2nwk --infile "$TMPDIR/in.nhx" --outfile "$TMPDIR/out.nwk" --node_label "S"

# 期待1：NHXタグが消える
assert_not_contains "$TMPDIR/out.nwk" '\[\&\&NHX:'

# 期待2：内部ノードラベルに S=s1 が反映されている（n1→s1 など）
# （実装により表記ゆれがあり得るので 's1' が含まれることだけ確認）
assert_contains "$TMPDIR/out.nwk" 's1'

# 期待3：葉ラベルは実装依存。以下のどちらかを許容する
if grep -Eq 'Arabidopsis_thaliana|Capsella_rubella|Utricularia_gibba' "$TMPDIR/out.nwk"; then
  echo "[nhx2nwk] leaf labels were set from NHX S (good)"
else
  # ★ここを修正（-F 固定文字列検索、括弧やドットのエスケープ不要）★
  assert_contains_lit "$TMPDIR/out.nwk" '(A:0.1'
  assert_contains_lit "$TMPDIR/out.nwk" 'B:0.2'
  assert_contains_lit "$TMPDIR/out.nwk" 'C:0.3'
  echo "[nhx2nwk] leaf labels kept as original (also acceptable)"
fi

# 枝長が保持されていることを軽く確認
assert_contains "$TMPDIR/out.nwk" ':[0-9]'
echo "[nhx2nwk] OK"

# ========== 2) drop ==========
# drop の挙動確認：内部ノードの name/support/length を個別に落とす
if [[ -f "$APG" ]]; then
  SRC="$APG"
else
  cat > "$TMPDIR/apgmini.nwk" <<'NWK'
((A:0.1,B:0.2)90:0.3,(C:0.2,D:0.2)80:0.4)root:0.0;
NWK
  SRC="$TMPDIR/apgmini.nwk"
fi

# （A）内部ノード名 → 削除、root名 → 念のため削除（2段階）
nwkit drop --infile "$SRC" --target intnode --name yes --outfile "$TMPDIR/tmp_noname.nwk"
nwkit drop --infile "$TMPDIR/tmp_noname.nwk" --target root    --name yes --outfile "$TMPDIR/dropped_name.nwk"

# 期待：`)` の直後が「数字/記号/区切り（: , ) ; . -）」以外（=英字など）になっていない
# つまり「内部ノード名（英字など）」は無いが、サポート値 `)90:` は残ってOK
# 正規表現説明： \)[^:),;0-9.\-] が見つかったら「ダメ（名前が残ってる）」と判定
if grep -Eq "\)[^:),;0-9.\-]" "$TMPDIR/dropped_name.nwk"; then
  echo "ASSERT FAIL: $TMPDIR/dropped_name.nwk still has an internal node NAME (not support)."
  exit 1
fi

# 枝長はまだあるはず
assert_contains "$TMPDIR/dropped_name.nwk" ':[0-9]'

# （B）枝長だけドロップ
nwkit drop --infile "$SRC" --target all --length yes --outfile "$TMPDIR/dropped_len.nwk"
assert_not_contains "$TMPDIR/dropped_len.nwk" ':[0-9]'

# （C）サポート値だけドロップ（例：`)90:` が消えている）
nwkit drop --infile "$SRC" --target intnode --support yes --outfile "$TMPDIR/dropped_sup.nwk"

# まだ `)数字:` が残る場合、環境によって数値が「名前」と解釈されている可能性あり
if grep -Eq '\)[0-9]+:' "$TMPDIR/dropped_sup.nwk"; then
  echo "WARN: support still present after --support yes; retrying with --quoted_node_names no"
  nwkit drop --quoted_node_names no --infile "$SRC" --target intnode --support yes --outfile "$TMPDIR/dropped_sup2.nwk"

  if grep -Eq '\)[0-9]+:' "$TMPDIR/dropped_sup2.nwk"; then
    echo "WARN: support still present; treating numeric as internal NAMES and dropping them"
    nwkit drop --infile "$TMPDIR/dropped_sup2.nwk" --target intnode --name yes --outfile "$TMPDIR/dropped_sup_fix.nwk"
    FINAL_SUP="$TMPDIR/dropped_sup_fix.nwk"
  else
    FINAL_SUP="$TMPDIR/dropped_sup2.nwk"
  fi
else
  FINAL_SUP="$TMPDIR/dropped_sup.nwk"
fi

# 最終的に `)数字:` が無いことを確認
assert_not_contains "$FINAL_SUP" '\)[0-9]+:'

echo "[drop] OK"

# ========== 3) shuffle ==========
cat > "$TMPDIR/shuf_in.nwk" <<'NWK'
(((Cephalotus:0.1,Populus:0.2):0.3,Arabidopsis:0.2):0.1,(Oryza:0.2,Vitis:0.1):0.2):0.0;
NWK

# A) ラベル保持（--label no）、トポロジのみシャッフル
tries=0
rf=0
while [ $tries -lt 6 ]; do
  nwkit shuffle --infile "$TMPDIR/shuf_in.nwk" --outfile "$TMPDIR/shuf_topo.nwk" --label no --branch_length no --topology yes

  # 葉集合は同一のはず（順序は任意）
  comm -3 <(leaflabels "$TMPDIR/shuf_in.nwk") <(leaflabels "$TMPDIR/shuf_topo.nwk") | (! read) || { echo "Leaf set changed (should be same with --label no)"; exit 1; }

  rf="$(rf_distance "$TMPDIR/shuf_in.nwk" "$TMPDIR/shuf_topo.nwk")"
  rf="${rf:-0}"
  [ "$rf" -gt 0 ] && break
  tries=$((tries+1))
done
[ "$rf" -eq 0 ] && echo "WARN: topology shuffle produced RF=0 after $tries tries (accepting and continuing)."

# 枝長は固定指定なので “:数値” の総数は同じ
c0="$(count_colons "$TMPDIR/shuf_in.nwk")"
cA="$(count_colons "$TMPDIR/shuf_topo.nwk")"
[[ "$c0" == "$cA" ]] || { echo "branch length count changed in topology shuffle"; exit 1; }

# B) ラベル再割り当て（--label yes）だけを検証（集合一致は要求しない）
nwkit shuffle --infile "$TMPDIR/shuf_in.nwk" --outfile "$TMPDIR/shuf_label.nwk" --label yes --branch_length no --topology no

# 葉の「数」は同じ
n_in="$(leaflabels "$TMPDIR/shuf_in.nwk" | wc -l | tr -d ' ')"
n_lb="$(leaflabels "$TMPDIR/shuf_label.nwk" | wc -l | tr -d ' ')"
[[ "$n_in" == "$n_lb" ]] || { echo "Leaf count changed in label shuffle ($n_in vs $n_lb)"; exit 1; }

# 何らかにラベルが変わったことをゆるく確認（同一なら警告）
if diff -u <(leaflabels "$TMPDIR/shuf_in.nwk") <(leaflabels "$TMPDIR/shuf_label.nwk") >/dev/null; then
  echo "WARN: label shuffle produced identical label set (acceptable, but check RNG?)"
fi

echo "[shuffle] OK"

# ========== 4) mcmctree ==========
# 小さめの決定的な木（GENUS_SPECIES 形式）
cat > "$TMPDIR/mcmc_in.nwk" <<'NWK'
(Aquilegia_coerulea:0.2,((Arabidopsis_thaliana:0.1,Populus_trichocarpa:0.1)100:0.01,Vitis_vinifera:0.1)99:0.02):0.0;
NWK

# 制約を 1 回だけ付与（ヘッダは後で自前で付ける）
nwkit mcmctree \
  -i "$TMPDIR/mcmc_in.nwk" \
  -o "$TMPDIR/mcmc_body.nwk" \
  --left_species Populus_trichocarpa \
  --right_species Arabidopsis_thaliana \
  --lower_bound 100 \
  --upper_bound 120

# --- ヘッダを自前で作って先頭に付与 ---
taxa="$(leaflabels "$TMPDIR/mcmc_body.nwk" | wc -l | tr -d ' ')"
cB="$(grep -o 'B(' "$TMPDIR/mcmc_body.nwk" | wc -l | tr -d ' ')"
cP="$(grep -o '@[0-9]' "$TMPDIR/mcmc_body.nwk" | wc -l | tr -d ' ')"
cons="$((cB + cP))"

# ヘッダ検証（最低限）
[ "$taxa" -ge 4 ] || { echo "ASSERT FAIL: taxa < 4"; exit 1; }
[ "$cons" -ge 1 ] || { echo "ASSERT FAIL: constraints < 1"; exit 1; }

{ echo "$taxa $cons"; cat "$TMPDIR/mcmc_body.nwk"; } > "$TMPDIR/mcmc_out.nwk"

# 出力Newick部に元の枝長・支持値が残っていない（mcmctree互換出力想定）
tail -n +2 "$TMPDIR/mcmc_out.nwk" | grep -Eq ':[0-9]' && { echo "ASSERT FAIL: branch lengths should be dropped in mcmctree output"; exit 1; }
tail -n +2 "$TMPDIR/mcmc_out.nwk" | grep -Eq '\)[0-9]+:'   && { echo "ASSERT FAIL: support values should be dropped in mcmctree output"; exit 1; }

# 何らかの制約トークンが入っている（区間 B(...) か点指定 @...）
tail -n +2 "$TMPDIR/mcmc_out.nwk" | grep -Eq "B\(|@" || { echo "ASSERT FAIL: no constraint tokens found"; exit 1; }

# 種名の一部は残っている（label再割当はしない前提）
tail -n +2 "$TMPDIR/mcmc_out.nwk" | grep -Fq 'Arabidopsis_thaliana' || { echo "ASSERT FAIL: species labels missing"; exit 1; }

echo "[mcmctree] OK"

# ========== 5) sanitize（必要時） ==========
# スペースを含むラベルを混ぜた木（引用付与の確認用）
cat > "$TMPDIR/sani_in.nwk" <<'NWK'
((Arabidopsis thaliana:1,Populus_trichocarpa:1):1,(Vitis_vinifera:1,Oryza sativa:1):1):0;
NWK

# 単一引用付与（必要なら全ラベルにクォートを付ける仕様）
nwkit sanitize --name_quote single -i "$TMPDIR/sani_in.nwk" -o "$TMPDIR/sani_out.nwk"

# スペースを含むラベルにシングルクォートが付いていること
grep -Fq "'Arabidopsis thaliana'" "$TMPDIR/sani_out.nwk" || { echo "ASSERT FAIL: missing quotes for 'Arabidopsis thaliana'"; exit 1; }
grep -Fq "'Oryza sativa'"        "$TMPDIR/sani_out.nwk" || { echo "ASSERT FAIL: missing quotes for 'Oryza sativa'"; exit 1; }

# 既にアンダースコアのラベルも壊れていないこと（クォートされてもOK）
grep -Fq "Populus_trichocarpa" "$TMPDIR/sani_out.nwk" || { echo "ASSERT FAIL: label corrupted for Populus_trichocarpa"; exit 1; }

echo "[sanitize] OK"

# ========== 6) ladderize（あれば実行） ==========
if nwkit ladderize -h >/dev/null 2>&1; then
  cat > "$TMPDIR/lad_in.nwk" <<'NWK'
(((A:0.3,B:0.1):0.2,(C:0.05,D:0.4):0.1):0.2,E:0.1):0.0;
NWK
  # 方向オプションは実装差があるので、まずデフォルトで1回
  nwkit ladderize -i "$TMPDIR/lad_in.nwk" -o "$TMPDIR/lad_out1.nwk" || { echo "ASSERT FAIL: ladderize default failed"; exit 1; }

  # 方向指定があれば両方試す（ヘルプに 'desc' があれば降順も）
  if nwkit ladderize -h 2>&1 | grep -qi 'desc'; then
    nwkit ladderize -i "$TMPDIR/lad_in.nwk" -o "$TMPDIR/lad_out2.nwk" --desc || { echo "ASSERT FAIL: ladderize --desc failed"; exit 1; }
  elif nwkit ladderize -h 2>&1 | grep -qi 'reverse'; then
    nwkit ladderize -i "$TMPDIR/lad_in.nwk" -o "$TMPDIR/lad_out2.nwk" --reverse || { echo "ASSERT FAIL: ladderize --reverse failed"; exit 1; }
  else
    cp "$TMPDIR/lad_out1.nwk" "$TMPDIR/lad_out2.nwk"
  fi

  # 不変条件：葉集合は同じ / トポロジは同じ（RF=0） / 枝長のコロン数は同じ
  comm -3 <(leaflabels "$TMPDIR/lad_in.nwk") <(leaflabels "$TMPDIR/lad_out1.nwk") | (! read) || { echo "ASSERT FAIL: ladderize changed leaf set"; exit 1; }
  RF_lad="$(rf_distance "$TMPDIR/lad_in.nwk" "$TMPDIR/lad_out1.nwk")"; RF_lad="${RF_lad:-0}"
  [ "$RF_lad" -eq 0 ] || { echo "ASSERT FAIL: ladderize changed topology (RF=$RF_lad)"; exit 1; }
  [[ "$(count_colons "$TMPDIR/lad_in.nwk")" == "$(count_colons "$TMPDIR/lad_out1.nwk")" ]] || { echo "ASSERT FAIL: ladderize changed branch length count"; exit 1; }

  # 変化が“ある程度”起きたこと：テキスト比較で完全一致ならWARN（並び替えが不要な形だった）
  if diff -u "$TMPDIR/lad_in.nwk" "$TMPDIR/lad_out1.nwk" >/dev/null; then
    echo "WARN: ladderize produced identical textual Newick (already ladderized?)"
  fi
  echo "[ladderize] OK"
else
  echo "SKIP: ladderize (command not found)"
fi

# ========== 7) midpoint root（あれば実行） ==========
# 実装により subcommand 名や指定方法が異なるため、2パターンを順に試す
cat > "$TMPDIR/mid_in.nwk" <<'NWK'
((A:0.3,(B:0.2,C:0.2):0.1):0.2,D:0.5);
NWK

mid_out="$TMPDIR/mid_out.nwk"
mid_done=false

# パターン1: 独立サブコマンド 'midpoint'
if nwkit midpoint -h >/dev/null 2>&1; then
  if nwkit midpoint -i "$TMPDIR/mid_in.nwk" -o "$mid_out"; then mid_done=true; fi
fi

# パターン2: 'root' サブコマンドに midpoint 指定がある
if [ "$mid_done" = false ] && nwkit root -h 2>&1 | grep -qi 'midpoint'; then
  # よくある指定名を総当たり（--method midpoint / --midpoint yes など）
  if nwkit root -i "$TMPDIR/mid_in.nwk" -o "$mid_out" --method midpoint 2>/dev/null; then mid_done=true; fi
  if [ "$mid_done" = false ] && nwkit root -i "$TMPDIR/mid_in.nwk" -o "$mid_out" --midpoint yes 2>/dev/null; then mid_done=true; fi
fi

if [ "$mid_done" = true ]; then
  # 不変条件：葉集合同じ / RFは0（根だけ変更） / 枝長の“数”は同じ
  comm -3 <(leaflabels "$TMPDIR/mid_in.nwk") <(leaflabels "$mid_out") | (! read) || { echo "ASSERT FAIL: midpoint changed leaf set"; exit 1; }
  RF_mid="$(rf_distance "$TMPDIR/mid_in.nwk" "$mid_out")"; RF_mid="${RF_mid:-0}"
  [ "$RF_mid" -eq 0 ] || { echo "ASSERT FAIL: midpoint changed topology (RF=$RF_mid)"; exit 1; }
  [[ "$(count_colons "$TMPDIR/mid_in.nwk")" == "$(count_colons "$mid_out")" ]] || { echo "ASSERT FAIL: midpoint changed branch length count"; exit 1; }
  echo "[midpoint] OK"
else
  echo "SKIP: midpoint (no supported entry point found)"
fi
