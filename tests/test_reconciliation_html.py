"""Exercise the shipped offline JavaScript with a deterministic DOM boundary."""

import copy
import json
import os
import shutil
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from nwkit.reconciliation_composition import comparison_scene
from nwkit.reconciliation_document import read_document, validate_document
from nwkit.reconciliation_html import SCRIPT, _svg
from nwkit.reconciliation_scene import build_scene

TOY = Path(__file__).parents[1] / "examples/reconciliation/toy.json"


def run_script(scene, source, view, actions):
    node = os.environ.get("NWKIT_TEST_NODE") or shutil.which("node")
    if node is None:
        pytest.skip("Node.js is unavailable for offline JavaScript behavior checks")
    svg = ET.fromstring(_svg(scene, 8))
    details = [
        json.loads(g[0].attrib["data-detail"])
        for g in svg
        if g.attrib.get("class") == "item"
    ]
    code = r"""
const assert=require('node:assert/strict'),vm=require('node:vm');
class Element {
 constructor(){this.attrs={};this.dataset={};this.listeners={};this.disabled=false;this.textContent='';this.classList={add(){},remove(){},toggle(){}};}
 addEventListener(k,fn){this.listeners[k]=fn;}
 setAttribute(k,v){this.attrs[k]=v;} removeAttribute(k){delete this.attrs[k];}
 getAttribute(k){return this.attrs[k];}
}
const elements=new Map();const get=id=>{if(!elements.has(id))elements.set(id,new Element());return elements.get(id);};
const details=DETAILS;
const entries=details.map(d=>{const e=new Element();e.firstElementChild=new Element();e.firstElementChild.dataset.detail=JSON.stringify(d);return e;});
get('source').textContent=JSON.stringify(SOURCE);get('viewstate').textContent=JSON.stringify(VIEW);
get('scene').attrs.viewBox='0 0 720 640';get('scene').dataset.full='0 0 720 900';get('scene').querySelectorAll=()=>entries;
const context={document:{getElementById:get},WeakMap,Map,Set,JSON,Uint8Array,console,setTimeout,assert};
vm.createContext(context);vm.runInContext(SCRIPT,context);vm.runInContext(ACTIONS,context);
console.log('Offline JavaScript checks passed');
"""
    replacements = {
        "DETAILS": details,
        "SOURCE": source,
        "VIEW": view,
        "SCRIPT": SCRIPT,
        "ACTIONS": actions,
    }
    for key, value in replacements.items():
        code = code.replace(key, json.dumps(value))
    try:
        result = subprocess.run([node], input=code, text=True, capture_output=True)
    except OSError as error:
        if os.environ.get("NWKIT_TEST_NODE"):
            raise
        pytest.skip(f"Optional Node.js cannot execute on this host: {error}")
    assert result.returncode == 0, result.stderr


def test_offline_collapse_expands_and_reset_restores_initial_geometry():
    doc = read_document(TOY)
    tree = doc.trees["lower:0"]
    view = {
        "schema_version": "0.1.0",
        "kind": "reconciliation_view",
        "layout": "paired",
        "collapsed_clades": [{"tree": tree.id, "node": tree.root}],
    }
    actions = f"""
assert(items.some(x=>Object.hasOwn(x.attrs,'hidden')));
const root=items.find(x=>details(x).tree==='lower:0'&&details(x).owner[1]==={json.dumps(tree.root)});
select(root);assert.equal(document.getElementById('collapse').disabled,false);
document.getElementById('collapse').onclick();assert(items.every(x=>!Object.hasOwn(x.attrs,'hidden')));
document.getElementById('collapse').onclick();assert(items.some(x=>Object.hasOwn(x.attrs,'hidden')));
document.getElementById('reset').onclick();assert(items.every(x=>!Object.hasOwn(x.attrs,'hidden')));
"""
    run_script(build_scene(doc, view=view), doc.record, view, actions)


def test_offline_details_are_scoped_to_panel_and_upper_collapse_stays_valid():
    first = read_document(TOY)
    second = copy.deepcopy(first.record)
    second["reconciliations"][0]["model"]["audit_panel"] = "second"
    docs = {"a": first, "b": validate_document(second)}
    runtime = {"seconds": None, "measured_at": None, "conditions": {}}
    view = {
        "schema_version": "0.1.0",
        "kind": "reconciliation_view",
        "layout": "topology_embedded",
    }
    record = {
        "kind": "reconciliation_comparison",
        "mode": "general",
        "view": view,
        "panels": [
            {
                "id": key,
                "label": key,
                "history_status": "available",
                "inference_status": "success",
                "document": key + ".json",
                "inference_runtime": runtime,
                "history_runtime": runtime,
            }
            for key in docs
        ],
    }
    source = dict(record, documents=[d.record for d in docs.values()])
    scene = comparison_scene(
        record, docs, width=720, font_size=8, apply_display_view=False
    )
    run_script(
        scene,
        source,
        view,
        r"""
select(items.find(x=>details(x).panel==='b'&&details(x).kind==='event'));
const records=JSON.parse(document.getElementById('detail').textContent).records;
assert(records.length>0);assert(records.every(r=>r.panel==='b'&&r.model.audit_panel==='second'));
select(items.find(x=>details(x).panel==='b'&&details(x).tree==='upper'));
assert.equal(document.getElementById('collapse').disabled,true);
document.getElementById('collapse').onclick();assert.equal(state.collapsed_clades,undefined);
""",
    )
