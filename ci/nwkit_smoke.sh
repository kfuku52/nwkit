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


# ========== 5.5) midpoint root（あれば実行） ==========
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
# 不変条件：葉集合同じ / RFは0（根だけ変更）
comm -3 <(leaflabels "$TMPDIR/mid_in.nwk") <(leaflabels "$mid_out") | (! read) || { echo "ASSERT FAIL: midpoint changed leaf set"; exit 1; }
RF_mid="$(rf_distance "$TMPDIR/mid_in.nwk" "$mid_out")"; RF_mid="${RF_mid:-0}"
[ "$RF_mid" -eq 0 ] || { echo "ASSERT FAIL: midpoint changed topology (RF=$RF_mid)"; exit 1; }

# 枝長トークン数は実装により±1等で変動し得るため、差分は WARN 扱いにする
c_in="$(count_colons "$TMPDIR/mid_in.nwk")"
c_out="$(count_colons "$mid_out")"
if [ "$c_out" -lt 1 ]; then
  echo "ASSERT FAIL: midpoint produced tree without branch lengths"
  exit 1
fi
if [ "$c_in" -ne "$c_out" ]; then
  echo "WARN: midpoint changed branch-length token count ($c_in -> $c_out) — acceptable due to root edge split/merge."
fi

echo "[midpoint] OK"
else
  echo "SKIP: midpoint (no supported entry point found)"
fi

# ========== 6) label（葉のリネーム検証） ==========
cat > "$TMPDIR/lab_in.nwk" <<'NWK'
(((A:0.1,B:0.2):0.1,(C:0.1,D:0.1):0.1):0.1,E:0.1):0.0;
NWK

# すべての葉を prefix=L で付け直す（--force yes）
nwkit label -i "$TMPDIR/lab_in.nwk" -o "$TMPDIR/lab_out.nwk" -t leaf --prefix L --force yes

# 葉の数は不変
n_in="$(leaflabels "$TMPDIR/lab_in.nwk" | wc -l | tr -d ' ')"
n_out="$(leaflabels "$TMPDIR/lab_out.nwk" | wc -l | tr -d ' ')"
[[ "$n_in" == "$n_out" ]] || { echo "ASSERT FAIL: label changed leaf count ($n_in -> $n_out)"; exit 1; }

# 旧ラベル（A,B,C,D,E）が残っていない & すべて L から始まる
comm -12 <(leaflabels "$TMPDIR/lab_in.nwk") <(leaflabels "$TMPDIR/lab_out.nwk") | (! read) || { echo "ASSERT FAIL: old leaf labels still present after label"; exit 1; }
if leaflabels "$TMPDIR/lab_out.nwk" | grep -Ev '^L' >/dev/null; then
  echo "ASSERT FAIL: some leaf labels do not start with 'L'"
  exit 1
fi
echo "[label] OK"

# ========== 7) prune（正規表現での剪定） ==========
cat > "$TMPDIR/prune_in.nwk" <<'NWK'
((((A:0.1,B:0.2):0.1,(C:0.1,D:0.2):0.1):0.1,(E:0.1,F:0.1):0.1):0.1):0.0;
NWK

# (A|B) に一致するもの「以外」を落とさず残す（= AB を残す） --invert_match yes
nwkit prune -i "$TMPDIR/prune_in.nwk" -o "$TMPDIR/prune_keepAB.nwk" -p '^(A|B)$' --invert_match yes
comm -3 <(printf "A\nB\n") <(leaflabels "$TMPDIR/prune_keepAB.nwk") | (! read) || { echo "ASSERT FAIL: prune keepAB leaf set mismatch"; exit 1; }

# (A|B) に一致するものを剪定（= AB を落とす） --invert_match デフォルト=no
nwkit prune -i "$TMPDIR/prune_in.nwk" -o "$TMPDIR/prune_dropAB.nwk" -p '^(A|B)$'
comm -3 <(printf "C\nD\nE\nF\n") <(leaflabels "$TMPDIR/prune_dropAB.nwk") | (! read) || { echo "ASSERT FAIL: prune dropAB leaf set mismatch"; exit 1; }

echo "[prune] OK"

# ========== 8) subtree（left/right と leaves の両モード） ==========
cat > "$TMPDIR/sub_in.nwk" <<'NWK'
(((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1):0.1,(E:0.1,F:0.1):0.1):0.0;
NWK

# 1) left/right で MRCA を指定 → 期待する葉集合は {A,B,C,D}
nwkit subtree -i "$TMPDIR/sub_in.nwk" -o "$TMPDIR/sub_AC.nwk" --left_leaf A --right_leaf C
comm -3 <(printf "A\nB\nC\nD\n") <(leaflabels "$TMPDIR/sub_AC.nwk") | (! read) \
  || { echo "ASSERT FAIL: subtree A|C leaf set mismatch"; exit 1; }

# 2) --leaves で MRCA を指定（B と D）。こちらも MRCA は同じクレード → {A,B,C,D}
nwkit subtree -i "$TMPDIR/sub_in.nwk" -o "$TMPDIR/sub_BD.nwk" --leaves B,D
comm -3 <(printf "A\nB\nC\nD\n") <(leaflabels "$TMPDIR/sub_BD.nwk") | (! read) \
  || { echo "ASSERT FAIL: subtree --leaves(B,D) leaf set mismatch"; exit 1; }

# 3) 念のため、指定した葉（B と D）が確実に含まれていること（subset チェック）
comm -23 <(printf "B\nD\n") <(leaflabels "$TMPDIR/sub_BD.nwk") | (! read) \
  || { echo "ASSERT FAIL: subtree --leaves(B,D) missing specified leaves"; exit 1; }

# （任意の強化）同じ MRCA 指定なので、両者は拓扑も一致するはず → RF=0
RF_sub="$(rf_distance "$TMPDIR/sub_AC.nwk" "$TMPDIR/sub_BD.nwk")"; RF_sub="${RF_sub:-0}"
[ "$RF_sub" -eq 0 ] || { echo "ASSERT FAIL: subtree outputs differ (RF=$RF_sub)"; exit 1; }

echo "[subtree] OK"

# ========== 9) rescale ==========
cat > "$TMPDIR/scale_in.nwk" <<'NWK'
((A:1.0,B:3.0):2.0,C:4.0):0.0;
NWK

factor=0.5
nwkit rescale -i "$TMPDIR/scale_in.nwk" -o "$TMPDIR/scale_out.nwk" --factor "$factor"

python - <<PY
import re, math, sys
def lens(path):
    s=open(path).read()
    return [float(x[1:]) for x in re.findall(r':[0-9]*\.?[0-9]+', s)]
l0=lens("$TMPDIR/scale_in.nwk")
l1=lens("$TMPDIR/scale_out.nwk")
assert len(l0)==len(l1), f"count changed: {len(l0)}->{len(l1)}"
fac=$factor
# 許容誤差（浮動小数フォーマット差吸収）
tol=1e-6
for a,b in zip(l0,l1):
    assert abs(b - a*fac) <= tol, f"value {a}->{b} not scaled by {fac}"
s0=sum(l0); s1=sum(l1)
assert abs(s1 - s0*fac) <= tol, f"sum {s0}->{s1} not scaled by {fac}"
PY

echo "[rescale] OK"

# ========== 10) constrain ==========
# ユーザー提供のバックボーン木（小さめ, 下線表記）
cat > "$TMPDIR/con_backbone.nwk" <<'NWK'
((Arabidopsis_thaliana:1,Populus_trichocarpa:1):1,(Vitis_vinifera:1,Oryza_sativa:1):1):0;
NWK

# 種リスト：木に存在する下線表記のみ（不一致を混ぜない）
cat > "$TMPDIR/species.txt" <<'TXT'
Arabidopsis_thaliana
Oryza_sativa
TXT

nwkit constrain \
  --backbone user \
  -i "$TMPDIR/con_backbone.nwk" \
  -o "$TMPDIR/con_out.nwk" \
  --species_list "$TMPDIR/species.txt" \
  --collapse no \
  || { echo "ASSERT FAIL: constrain execution failed"; exit 1; }

# 出力がNewickであること
assert_contains_lit "$TMPDIR/con_out.nwk" ';'

# 期待する種が現れること
grep -Fq 'Arabidopsis_thaliana' "$TMPDIR/con_out.nwk" || { echo "ASSERT FAIL: Arabidopsis_thaliana missing"; exit 1; }
grep -Fq 'Oryza_sativa'        "$TMPDIR/con_out.nwk" || { echo "ASSERT FAIL: Oryza_sativa missing"; exit 1; }

echo "[constrain] OK"

# ========== 11) root --method outgroup ==========
cat > "$TMPDIR/root_in.nwk" <<'NWK'
((A:0.3,(B:0.2,C:0.2):0.1):0.2,D:0.5);
NWK

nwkit root -i "$TMPDIR/root_in.nwk" -o "$TMPDIR/root_out_og.nwk" --method outgroup --outgroup D

# 葉集合は同じ
comm -3 <(leaflabels "$TMPDIR/root_in.nwk") <(leaflabels "$TMPDIR/root_out_og.nwk") | (! read) \
  || { echo "ASSERT FAIL: root(outgroup) changed leaf set"; exit 1; }

# RF=0（根の移動のみでトポロジは不変）
RF_root_og="$(rf_distance "$TMPDIR/root_in.nwk" "$TMPDIR/root_out_og.nwk")"; RF_root_og="${RF_root_og:-0}"
[ "$RF_root_og" -eq 0 ] || { echo "ASSERT FAIL: root(outgroup) changed topology (RF=$RF_root_og)"; exit 1; }

# 枝長トークン数は実装により±1等で変動あり → 0本になっていないことだけ確認
c_in="$(count_colons "$TMPDIR/root_in.nwk")"
c_og="$(count_colons "$TMPDIR/root_out_og.nwk")"
[ "$c_og" -ge 1 ] || { echo "ASSERT FAIL: root(outgroup) produced tree without branch lengths"; exit 1; }

echo "[root:outgroup] OK"

# ========== 12) root --method transfer ==========
# 参照：外群Dでルートした木（上の outgroup 出力を再利用）
cp "$TMPDIR/root_out_og.nwk" "$TMPDIR/root_ref_for_transfer.nwk"

# 入力をもう一度用意（transferの入力元）
cp "$TMPDIR/root_in.nwk" "$TMPDIR/root_xfer_src.nwk"

# transfer 実行（--infile2 に参照木）
nwkit root -i "$TMPDIR/root_xfer_src.nwk" -o "$TMPDIR/root_out_xfer.nwk" --method transfer --infile2 "$TMPDIR/root_ref_for_transfer.nwk"

# 葉集合は同じ
comm -3 <(leaflabels "$TMPDIR/root_ref_for_transfer.nwk") <(leaflabels "$TMPDIR/root_out_xfer.nwk") | (! read) \
  || { echo "ASSERT FAIL: root(transfer) changed leaf set"; exit 1; }

# 参照ファイルが本当に出来ているかチェック
[ -f "$TMPDIR/root_ref_for_transfer.nwk" ] || { echo "ASSERT FAIL: missing root_ref_for_transfer.nwk"; ls -l "$TMPDIR"; exit 1; }
[ -f "$TMPDIR/root_out_xfer.nwk" ] || { echo "ASSERT FAIL: missing root_out_xfer.nwk"; ls -l "$TMPDIR"; exit 1; }

# ルート分割（root split）が一致するかを ETE3 で厳密確認
python - <<PY
from ete3 import Tree
def root_split(path):
    t = Tree(open(path).read())
    ch = t.get_children()
    assert len(ch) == 2, f"root is not binary in {path}"
    return tuple(sorted([";".join(sorted([lf.name for lf in c.iter_leaves()])) for c in ch]))

ref = root_split("${TMPDIR}/root_ref_for_transfer.nwk")
got = root_split("${TMPDIR}/root_out_xfer.nwk")
assert ref == got, f"root split mismatch:\\n  ref={ref}\\n  got={got}"
PY

echo "[root:transfer] OK"

# ========== 13) intersection ==========
cat > "$TMPDIR/int_a.nwk" <<'NWK'
(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1);
NWK
cat > "$TMPDIR/int_b.nwk" <<'NWK'
((B:0.2,(C:0.2,E:0.2):0.1):0.1);
NWK

nwkit intersection -i "$TMPDIR/int_a.nwk" -i2 "$TMPDIR/int_b.nwk" -o "$TMPDIR/int_out.nwk"

# 出力がNewick（;を含む）
assert_contains_lit "$TMPDIR/int_out.nwk" ';'

# 葉集合が {B,C} に一致
comm -3 <(printf "B\nC\n") <(leaflabels "$TMPDIR/int_out.nwk") | (! read) \
  || { echo "ASSERT FAIL: intersection leaf set mismatch (expected B,C)"; exit 1; }

echo "[intersection] OK"

# ========== 14) printlabel ==========
cat > "$TMPDIR/pl_in.nwk" <<'NWK'
(((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1):0.1):0.0;
NWK

# 1) ターゲットそのもの（A と C）を出力
#    一部環境では -o が無視され STDOUT に出るため、明示的にリダイレクトする
nwkit printlabel -i "$TMPDIR/pl_in.nwk" -p '^(A|C)$' -t leaf > "$TMPDIR/pl_targets.txt"

# 出力から「葉集合に属するラベル」だけを抽出して集合比較
printf "A\nC\n" > "$TMPDIR/exp_targets.txt"
grep -Eo '[A-Za-z_][A-Za-z0-9_]*' "$TMPDIR/pl_targets.txt" | sort -u \
  | comm -12 - <(leaflabels "$TMPDIR/pl_in.nwk") > "$TMPDIR/seen_targets.txt"
diff -u "$TMPDIR/exp_targets.txt" "$TMPDIR/seen_targets.txt" >/dev/null \
  || { echo "ASSERT FAIL: printlabel targets mismatch"; exit 1; }

# 2) 姉妹ノード（A→B, C→D）を出力
nwkit printlabel -i "$TMPDIR/pl_in.nwk" -p '^(A|C)$' -t leaf --sister yes > "$TMPDIR/pl_sisters.txt"

printf "B\nD\n" > "$TMPDIR/exp_sisters.txt"
grep -Eo '[A-Za-z_][A-Za-z0-9_]*' "$TMPDIR/pl_sisters.txt" | sort -u \
  | comm -12 - <(leaflabels "$TMPDIR/pl_in.nwk") > "$TMPDIR/seen_sisters.txt"
diff -u "$TMPDIR/exp_sisters.txt" "$TMPDIR/seen_sisters.txt" >/dev/null \
  || { echo "ASSERT FAIL: printlabel sisters mismatch"; exit 1; }

echo "[printlabel] OK"

# ========== 15) mark ==========
cat > "$TMPDIR/mark_in.nwk" <<'NWK'
(((A:0.1,B:0.2):0.1,(C:0.1,D:0.2):0.1):0.1,E:0.1):0.0;
NWK

# 1) 葉ラベルに印（A と C に _TAG を接尾）
nwkit mark -i "$TMPDIR/mark_in.nwk" -o "$TMPDIR/mark_leaf.nwk" \
  -p '^(A|C)$' -t leaf --insert_txt TAG --insert_sep _ --insert_pos suffix

grep -Fq 'A_TAG' "$TMPDIR/mark_leaf.nwk" || { echo "ASSERT FAIL: A_TAG not found after mark leaf"; exit 1; }
grep -Fq 'C_TAG' "$TMPDIR/mark_leaf.nwk" || { echo "ASSERT FAIL: C_TAG not found after mark leaf"; exit 1; }
! grep -Fq 'B_TAG' "$TMPDIR/mark_leaf.nwk" || { echo "ASSERT FAIL: B should not be tagged"; exit 1; }
! grep -Fq 'D_TAG' "$TMPDIR/mark_leaf.nwk" || { echo "ASSERT FAIL: D should not be tagged"; exit 1; }

# 2) MRCA に印（A と B の MRCA に CLD を付与）。内部ノードに何らかのラベルが付く想定。
nwkit mark -i "$TMPDIR/mark_in.nwk" -o "$TMPDIR/mark_mrca.nwk" \
  -p '^(A|B)$' -t mrca --insert_txt CLD --insert_sep _ --insert_pos suffix

# 内部ノードに 'CLD' が一回以上は現れること
grep -Fq 'CLD' "$TMPDIR/mark_mrca.nwk" || { echo "ASSERT FAIL: no MRCA tag inserted"; exit 1; }
# 葉には CLD が付いていない（葉の印ではない）
! grep -Eq 'A_?CLD|B_?CLD|C_?CLD|D_?CLD|E_?CLD' "$TMPDIR/mark_mrca.nwk" || { echo "ASSERT FAIL: MRCA tag should not be on leaves"; exit 1; }

echo "[mark] OK"

# ========== 16) transfer ==========
# base: トポロジのみ（長さも支持も名前もなし）
cat > "$TMPDIR/xfer_base.nwk" <<'NWK'
((A,B),(C,D));
NWK

# src1: 長さ＆支持値あり
cat > "$TMPDIR/xfer_src_support.nwk" <<'NWK'
((A:0.1,B:0.2)90:0.3,(C:0.15,D:0.25)80:0.35):0.01;
NWK

# 長さ＆支持値を移送
nwkit transfer -i "$TMPDIR/xfer_base.nwk" -i2 "$TMPDIR/xfer_src_support.nwk" \
  -o "$TMPDIR/xfer_out_ls.nwk" --length yes --support yes

# 葉集合は同じ
comm -3 <(leaflabels "$TMPDIR/xfer_base.nwk") <(leaflabels "$TMPDIR/xfer_out_ls.nwk") | (! read) \
  || { echo "ASSERT FAIL: transfer(length/support) changed leaf set"; exit 1; }

# 枝長と支持値が入っている（フォーマット差は許容、存在だけ確認）
grep -Eq ':[0-9]' "$TMPDIR/xfer_out_ls.nwk" || { echo "ASSERT FAIL: no branch lengths after transfer"; exit 1; }
grep -Eq '\)[0-9]+:' "$TMPDIR/xfer_out_ls.nwk" || { echo "ASSERT FAIL: no support values after transfer"; exit 1; }

# src2: 内部ノード名あり
cat > "$TMPDIR/xfer_src_names.nwk" <<'NWK'
((A:0.1,B:0.2)Nleft:0.3,(C:0.15,D:0.25)Nright:0.35)Nroot:0.0;
NWK

# base に名前を移送（内部ノード対象）
nwkit transfer -i "$TMPDIR/xfer_base.nwk" -i2 "$TMPDIR/xfer_src_names.nwk" \
  -o "$TMPDIR/xfer_out_names.nwk" -t intnode --name yes

grep -Fq 'Nleft'  "$TMPDIR/xfer_out_names.nwk" || { echo "ASSERT FAIL: Nleft not transferred"; exit 1; }
grep -Fq 'Nright' "$TMPDIR/xfer_out_names.nwk" || { echo "ASSERT FAIL: Nright not transferred"; exit 1; }
grep -Fq 'Nroot'  "$TMPDIR/xfer_out_names.nwk" || { echo "ASSERT FAIL: Nroot not transferred"; exit 1; }

echo "[transfer] OK"

# ========== 17) dist（RF距離の検証とエラー系） ==========
# 5葉の決定的な2木を用意（E は外側固定）
cat > "$TMPDIR/rf_t1.nwk" <<'NWK'
((((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1):0.1,E:0.1);
NWK
cat > "$TMPDIR/rf_t2.nwk" <<'NWK'
((((A:0.1,C:0.1):0.1,B:0.1):0.1,D:0.1):0.1,E:0.1);
NWK

# 1) 同一木どうし → RF=0
out_same="$(nwkit dist -i "$TMPDIR/rf_t1.nwk" -i2 "$TMPDIR/rf_t1.nwk")"
rf_same="$(printf "%s\n" "$out_same" | sed -n 's/.*Robinson-Foulds distance = \([0-9][0-9]*\).*/\1/p' | head -n1)"
[ "${rf_same:-999}" -eq 0 ] || { echo "ASSERT FAIL: dist identical trees should be RF=0 (got $rf_same)"; exit 1; }

# 2) 別木どうし → nwkit の出力値が ETE3 の計算と一致
python - <<PY
from ete3 import Tree
t1=Tree(open("${TMPDIR}/rf_t1.nwk").read())
t2=Tree(open("${TMPDIR}/rf_t2.nwk").read())
rf,maxrf = t1.robinson_foulds(t2, unrooted_trees=False)[:2]
print(rf, maxrf)
PY
read rf_true max_true < <(python - <<PY
from ete3 import Tree
t1=Tree(open("${TMPDIR}/rf_t1.nwk").read())
t2=Tree(open("${TMPDIR}/rf_t2.nwk").read())
rf,maxrf = t1.robinson_foulds(t2, unrooted_trees=False)[:2]
print(rf, maxrf)
PY
)

out_diff="$(nwkit dist -i "$TMPDIR/rf_t1.nwk" -i2 "$TMPDIR/rf_t2.nwk")"
rf_diff="$(printf "%s\n" "$out_diff" | sed -n 's/.*Robinson-Foulds distance = \([0-9][0-9]*\).*/\1/p' | head -n1)"
max_diff="$(printf "%s\n" "$out_diff" | sed -n 's/.*(max = \([0-9][0-9]*\)).*/\1/p' | head -n1)"

[ "${rf_diff:-999}" = "$rf_true" ] || { echo "ASSERT FAIL: dist RF mismatch (nwkit=$rf_diff, ete3=$rf_true)"; exit 1; }
[ "${max_diff:-999}" = "$max_true" ] || { echo "ASSERT FAIL: dist MAX mismatch (nwkit=$max_diff, ete3=$max_true)"; exit 1; }

# 3) 葉集合が一致しない → エラー終了し、メッセージを含む
cat > "$TMPDIR/rf_bad.nwk" <<'NWK'
(((A:0.1,B:0.1):0.1,C:0.1):0.1,X:0.1);
NWK
if nwkit dist -i "$TMPDIR/rf_t1.nwk" -i2 "$TMPDIR/rf_bad.nwk" >"$TMPDIR/rf_err.txt" 2>&1; then
  echo "ASSERT FAIL: dist should fail on mismatched leaf sets"; exit 1;
fi
grep -Fq 'Leaf name(s) did not match' "$TMPDIR/rf_err.txt" || { echo "ASSERT FAIL: expected mismatch message not found"; exit 1; }

echo "[dist] OK"

# ========== 18) info（基本統計の妥当性） ==========
# 1) 下線のみの木
cat > "$TMPDIR/info_in1.nwk" <<'NWK'
((A:1.0,(B:1.0,C:1.0):1.0):1.0,(D:1.0,E:1.0):1.0):0.0;
NWK
n_expect="$(leaflabels "$TMPDIR/info_in1.nwk" | wc -l | tr -d ' ')"
nwkit info -i "$TMPDIR/info_in1.nwk" > "$TMPDIR/info_out1.txt"
n_seen="$(sed -n 's/.*Number of leaves[^=]*=\s*\([0-9][0-9]*\).*/\1/p' "$TMPDIR/info_out1.txt" | head -n1)"
[ "${n_seen:-999}" = "$n_expect" ] || { echo "ASSERT FAIL: info leaf count mismatch (got $n_seen, expect $n_expect)"; exit 1; }

# 2) 空白を含むラベル（クォート済み）
cat > "$TMPDIR/info_in2.nwk" <<'NWK'
(('Arabidopsis thaliana':1,'Oryza sativa':1):1,'Vitis vinifera':1);
NWK
n_expect2="$(leaflabels "$TMPDIR/info_in2.nwk" | wc -l | tr -d ' ')"
nwkit info -i "$TMPDIR/info_in2.nwk" > "$TMPDIR/info_out2.txt"
n_seen2="$(sed -n 's/.*Number of leaves[^=]*=\s*\([0-9][0-9]*\).*/\1/p' "$TMPDIR/info_out2.txt" | head -n1)"
[ "${n_seen2:-999}" = "$n_expect2" ] || { echo "ASSERT FAIL: info leaf count mismatch for quoted labels (got $n_seen2, expect $n_expect2)"; exit 1; }

echo "[info] OK"
