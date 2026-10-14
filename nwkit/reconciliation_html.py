"""Self-contained, offline SVG viewer with original vector PDF and display state."""

from __future__ import annotations

import base64
import copy
import html
import io
import json
import math
import tempfile
import textwrap
import zipfile
from pathlib import Path

from nwkit.reconciliation_render import write_scene


def _json(value):
    return (
        json.dumps(value, ensure_ascii=True, allow_nan=False)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


def _svg(scene, font_size, font_family="DejaVu Sans"):
    font_family = html.escape(font_family, quote=True)
    elements = []
    for p in scene.primitives:
        data = html.escape(
            _json(
                {
                    "id": p.id,
                    "owner": p.owner,
                    "panel": p.panel,
                    "tree": p.tree,
                    "kind": p.kind,
                    "text": p.text,
                }
            ),
            quote=True,
        )
        attrs = f'data-detail="{data}" stroke-width="{p.width}"'
        color = html.escape(p.color)
        coords = " ".join(f"{x:.5f},{y:.5f}" for x, y in p.points)
        if p.kind in {"line", "upper"}:
            element = (
                f'<polyline points="{coords}" {attrs} stroke="{color}" fill="none"/>'
            )
        elif p.kind == "transfer":
            (x, y), (xx, yy) = p.points[0], p.points[-1]
            cx, cy = (x + xx) / 2 - (yy - y) * 0.25, (y + yy) / 2 + (xx - x) * 0.25
            element = f'<path d="M{x},{y} Q{cx},{cy} {xx},{yy}" marker-end="url(#arrow)" {attrs} stroke="{color}" fill="none"/>'
        elif p.kind == "event":
            x, y = p.points[0]
            element = _event_symbol(p.marker, x, y, attrs, color)
        else:
            x, y = p.points[0]
            lines = textwrap.wrap(p.text, 24) or [""]
            spans = "".join(
                f'<tspan x="{x}" y="{y + (index - (len(lines) - 1) / 2) * font_size * 1.2}">{html.escape(line)}</tspan>'
                for index, line in enumerate(lines)
            )
            element = f'<text {attrs} stroke="none" fill="{color}" font-size="{font_size}" font-family="{font_family}" dominant-baseline="middle">{spans}</text>'

        element = element.replace(" ", ' class="visual" ', 1)
        hit = ""
        if p.kind in {"line", "upper", "transfer"}:
            hit = (
                element.replace('class="visual"', 'class="hit-area"')
                .replace(f'stroke="{color}"', 'stroke="transparent"')
                .replace(f'stroke-width="{p.width}"', 'stroke-width="12"')
                .replace('marker-end="url(#arrow)"', "")
            )
        elements.append(
            f'<g class="item" tabindex="0" role="button">{hit}{element}</g>'
        )
    return f'<svg id="scene" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {scene.width} {min(scene.height + 20, 640)}" data-full="0 0 {scene.width} {scene.height + 20}" role="img" aria-label="Recorded reconciliation trees"><defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/></marker></defs>{"".join(elements)}</svg>'


def _event_symbol(marker, x, y, attrs, color):
    style = f'{attrs} stroke="{color}" fill="{color}"'
    if marker == "s":
        return f'<rect x="{x - 3}" y="{y - 3}" width="6" height="6" {style}/>'
    if marker in {"x", "X"}:
        return f'<path d="M{x - 3},{y - 3} L{x + 3},{y + 3} M{x - 3},{y + 3} L{x + 3},{y - 3}" {attrs} stroke="{color}" fill="none"/>'
    if marker in {">", "<"}:
        direction = 1 if marker == ">" else -1
        points = f"{x + 3 * direction},{y} {x - 3 * direction},{y - 3} {x - 3 * direction},{y + 3}"
        return f'<polygon points="{points}" {style}/>'
    if marker == "D":
        return f'<polygon points="{x},{y - 4} {x + 4},{y} {x},{y + 4} {x - 4},{y}" {style}/>'
    if marker == "*":
        points = " ".join(
            f"{x + (4 if i % 2 == 0 else 1.8) * math.cos(i * math.pi / 5 - math.pi / 2)},{y + (4 if i % 2 == 0 else 1.8) * math.sin(i * math.pi / 5 - math.pi / 2)}"
            for i in range(10)
        )
        return f'<polygon points="{points}" {style}/>'
    return f'<circle cx="{x}" cy="{y}" r="3" {style}/>'


SCRIPT = r"""
const source=JSON.parse(document.getElementById('source').textContent);
const state=JSON.parse(document.getElementById('viewstate').textContent);
const svg=document.getElementById('scene'), items=[...svg.querySelectorAll('.item')];
const initial=svg.getAttribute('viewBox');const full=svg.dataset.full;let selected=null,drag=null;
function details(item){return JSON.parse(item.firstElementChild.dataset.detail);}
function download(text,name,type){const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([text],{type}));a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);}
function documents(){return source.kind==='reconciliation_document'?[source]:source.documents;}
function documentFor(info){if(source.kind==='reconciliation_document')return source;const panels=source.panels.filter(p=>p.document!==null);const index=panels.findIndex(p=>p.id===info.panel);return source.documents[index]||null;}
function treeFor(info){const doc=documentFor(info);if(!doc||!info.owner)return null;let tree=doc.trees.find(t=>t.id===(info.tree||info.owner[0]));if(!tree){const rec=doc.reconciliations.find(r=>r.id===info.owner[0]);if(rec)tree=doc.trees.find(t=>t.id===rec.lower_tree);}return tree&&tree.nodes.some(n=>n.id===info.owner[1])?tree:null;}
function nodeSet(tree,node){const result=new Set();const nodes=new Map(tree.nodes.map(n=>[n.id,n]));let queue=[node];while(queue.length){let id=queue.pop();result.add(id);queue.push(...nodes.get(id).children);}return result;}
const cladeCache=new WeakMap();
function cladeKey(tree,node){let keys=cladeCache.get(tree);if(!keys){keys=new Map();const nodes=new Map(tree.nodes.map(n=>[n.id,n]));const order=[],stack=[tree.root],leaves=new Map();while(stack.length){const id=stack.pop();order.push(id);stack.push(...nodes.get(id).children);}for(const id of order.reverse()){const n=nodes.get(id);const labels=n.children.length?n.children.flatMap(c=>leaves.get(c)).sort():[n.label];leaves.set(id,labels);keys.set(id,JSON.stringify([tree.id,labels]));}cladeCache.set(tree,keys);}return keys.get(node);}
function collapseProblem(tree,node){if(state.layout==='topology_embedded'&&documents().some(doc=>doc.reconciliations.some(r=>r.upper_tree===tree.id)))return 'Upper context collapse requires a standalone or paired view.';const keys=documents().flatMap(doc=>doc.trees.filter(t=>t.id===tree.id&&t.nodes.some(n=>n.id===node)).map(t=>cladeKey(t,node)));return new Set(keys).size>1?'This tree/node ID names different clades across panels; a shared collapse would be ambiguous.':null;}
function applyCollapse(){for(const item of items)item.removeAttribute("hidden");for(const clade of state.collapsed_clades||[]){for(const item of items){const info=details(item);const tree=treeFor(info);if(tree&&tree.id===clade.tree&&tree.nodes.some(n=>n.id===clade.node)&&info.owner[1]!==clade.node&&nodeSet(tree,clade.node).has(info.owner[1]))item.setAttribute("hidden","");}}}
function select(item){items.forEach(x=>x.classList.remove('selected'));selected=details(item);const tree=treeFor(selected);for(const x of items){const d=details(x),other=treeFor(d);if(tree&&other&&cladeKey(tree,selected.owner[1])===cladeKey(other,d.owner[1]))x.classList.add('selected');}const metadata=[],doc=documentFor(selected);if(tree&&doc){for(const rec of doc.reconciliations){const lower=rec.lower_tree===tree.id,upper=rec.upper_tree===tree.id;if(!lower&&!upper)continue;for(const event of rec.events)if((lower?event.lower_node:event.upper_node)===selected.owner[1])metadata.push({panel:selected.panel,reconciliation:rec.id,history_kind:rec.history_kind,selection:rec.selection,model:rec.model,event});}}const button=document.getElementById('collapse');const problem=tree?collapseProblem(tree,selected.owner[1]):'Select a recorded clade.';button.disabled=!!problem;button.title=problem||'';document.getElementById('detail').textContent=JSON.stringify({primitive:selected,records:metadata,collapse_notice:problem},null,2);}
items.forEach(item=>{item.addEventListener('click',()=>select(item));item.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();select(item);}});});
document.getElementById('search').addEventListener('input',event=>{const q=event.target.value.toLowerCase();items.forEach(item=>{const d=details(item);item.classList.toggle('match',!!q&&JSON.stringify(d).toLowerCase().includes(q));});});
document.getElementById('collapse').onclick=()=>{if(!selected)return;const tree=treeFor(selected);if(!tree||collapseProblem(tree,selected.owner[1]))return;state.collapsed_clades??=[];const key={tree:tree.id,node:selected.owner[1]};let i=state.collapsed_clades.findIndex(x=>x.tree===key.tree&&x.node===key.node);if(i<0)state.collapsed_clades.push(key);else state.collapsed_clades.splice(i,1);applyCollapse();};
function zoom(scale){let [x,y,w,h]=svg.getAttribute('viewBox').split(' ').map(Number);if(w*scale<20||w*scale>100000)return;svg.setAttribute('viewBox',`${x+w*(1-scale)/2} ${y+h*(1-scale)/2} ${w*scale} ${h*scale}`);}
document.getElementById('zoomin').onclick=()=>zoom(1/1.3);
document.getElementById('zoomout').onclick=()=>zoom(1.3);
document.getElementById('fitall').onclick=()=>svg.setAttribute('viewBox',full);
document.getElementById('reset').onclick=()=>{svg.setAttribute('viewBox',initial);state.collapsed_clades=[];applyCollapse();};
document.getElementById('saveview').onclick=()=>download(JSON.stringify(state,null,2),'reconciliation.view.json','application/json');
document.getElementById('savesource').onclick=()=>{if(source.kind==='reconciliation_document'){download(JSON.stringify(source,null,2),'reconciliation.source.json','application/json');return;}const bytes=Uint8Array.from(atob(document.getElementById('sourcebundle').textContent),c=>c.charCodeAt(0));download(bytes,'reconciliation.sources.zip','application/zip');};
document.getElementById('print').onclick=()=>window.print();
document.getElementById('pdf').onclick=()=>{const bytes=Uint8Array.from(atob(document.getElementById('pdfbytes').textContent),c=>c.charCodeAt(0));download(bytes,'reconciliation.pdf','application/pdf');};
svg.addEventListener('wheel',event=>{event.preventDefault();let [x,y,w,h]=svg.getAttribute('viewBox').split(' ').map(Number);const scale=event.deltaY>0?1.15:1/1.15;if(w*scale<20||w*scale>100000)return;svg.setAttribute('viewBox',`${x+w*(1-scale)/2} ${y+h*(1-scale)/2} ${w*scale} ${h*scale}`);},{passive:false});
svg.addEventListener('pointerdown',event=>{if(event.target===svg){drag={x:event.clientX,y:event.clientY,box:svg.getAttribute('viewBox').split(' ').map(Number)};svg.setPointerCapture(event.pointerId);}});
svg.addEventListener('pointermove',event=>{if(!drag)return;const [x,y,w,h]=drag.box,rect=svg.getBoundingClientRect();svg.setAttribute('viewBox',`${x-(event.clientX-drag.x)*w/rect.width} ${y-(event.clientY-drag.y)*h/rect.height} ${w} ${h}`);});
svg.addEventListener('pointerup',()=>drag=null);applyCollapse();
"""


def write_html(
    scene,
    path,
    *,
    source,
    view,
    title,
    font_size,
    font_family,
    pdf_pages,
    pdf_page_size,
    legend,
    transparent=False,
    pdf_scene=None,
):
    with tempfile.TemporaryDirectory(prefix="nwkit-offline-pdf-") as folder:
        pdf = Path(folder) / "document.pdf"
        write_scene(
            pdf_scene if pdf_scene is not None else scene,
            pdf,
            title=title,
            image_format="pdf",
            font_size=font_size,
            font_family=font_family,
            pdf_pages=pdf_pages,
            pdf_page_size=pdf_page_size,
            legend=legend,
            transparent=transparent,
        )
        pdfbytes = base64.b64encode(pdf.read_bytes()).decode()
    bundle = ""
    if source.get("kind") == "reconciliation_comparison":
        manifest = copy.deepcopy(source)
        manifest.pop("documents", None)
        raw = manifest.pop("document_bytes")
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            index = 0
            for panel in manifest["panels"]:
                if panel["document"] is None:
                    continue
                filename = f"panel-{index}.document{Path(panel['document']).suffix}"
                panel["document"] = filename
                archive.writestr(filename, base64.b64decode(raw[index]))
                index += 1
            archive.writestr(
                "comparison.json", json.dumps(manifest, indent=2, allow_nan=False)
            )
        bundle = base64.b64encode(buffer.getvalue()).decode()
    source_button = "Download source bundle" if bundle else "Save source JSON"
    notes = "".join(f"<li>{html.escape(note)}</li>" for note in scene.notices)
    title = html.escape(title or "Recorded reconciliation")
    content = f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>{title}</title>
<style>body{{font:14px system-ui;color:#27364b;margin:18px}}header{{display:flex;gap:8px;flex-wrap:wrap}}button,input{{padding:7px}}main{{display:grid;grid-template-columns:minmax(0,1fr) 320px;gap:16px}}svg{{width:100%;min-height:500px;max-height:80vh;border:1px solid #d7dfe8;touch-action:none}}pre{{white-space:pre-wrap;word-break:break-word;max-height:80vh;overflow:auto}}.item[hidden]{{display:none}}.selected>.visual{{stroke:#f59e0b!important;stroke-width:3!important}}.match>.visual{{stroke:#e11d48!important;stroke-width:2!important}}@media print{{header,aside{{display:none}}main{{display:block}}svg{{max-height:none;min-height:0;border:0;width:100%;height:auto}}@page{{size:landscape;margin:12mm}}}}</style>
<h1>{title}</h1><p>{html.escape(scene.caption)}</p><header><input id="search" aria-label="Search labels, IDs and event kinds" placeholder="Search labels / IDs"><button id="collapse" disabled>Toggle selected clade</button><button id="reset">Reset view</button><button id="zoomin">Zoom in</button><button id="zoomout">Zoom out</button><button id="fitall">Fit all</button><button id="saveview">Save view JSON</button><button id="savesource">{source_button}</button><button id="pdf">Download vector PDF</button><button id="print">Print current view / PDF</button></header>
<p>The downloadable PDF preserves the initial saved view and pagination. Print uses the current browser display. Drag the background to pan; scroll or use zoom controls. Source data, histories and fitted parameters stay unchanged.</p><main><section>{_svg(scene, font_size, font_family)}<ul>{notes}</ul></section><aside><h2>Recorded details</h2><pre id="detail">Select a branch or event.</pre></aside></main>
<script type="application/json" id="source">{_json(source)}</script><script type="application/json" id="viewstate">{_json(view)}</script><script type="application/octet-stream" id="sourcebundle">{bundle}</script><script type="application/octet-stream" id="pdfbytes">{pdfbytes}</script><script>{SCRIPT}</script></html>"""
    Path(path).write_text(content, encoding="utf-8")
