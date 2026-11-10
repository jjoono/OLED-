"""Build a self-contained interactive 3D viewer (3Dmol.js) of all optimized structures."""
import os, json, glob

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")

# system tag -> (display name, E_b eV or None, note)
E = {}
for f in ["psi4_binding_refined_eV.json", "bphen_cs_binding_eV.json",
          "newcand_binding_eV.json", "htl_binding_eV.json", "zns_binding_eV.json",
          "b3pympm_binding_eV.json", "psi4_binding_eV.json"]:
    p = os.path.join(RUNS, f)
    if os.path.exists(p):
        for k, v in json.load(open(p)).items():
            if isinstance(v, float) and not k.endswith("_offset") and k not in E:
                E[k] = v

SYSTEMS = [
    ("HATCN_Ag_CN",    "HATCN + Ag (nitrile N)",       E.get("HATCN_Ag_CN")),
    ("HATCN_Ag_face",  "HATCN + Ag (face)",            E.get("HATCN_Ag_face")),
    ("TPBi_Ag",        "TPBi + Ag (benzimidazole N)",  E.get("TPBi_Ag")),
    ("pbPPhenB_Ag",    "p-bPPhenB + Ag (chelate)",     E.get("pbPPhenB_Ag")),
    ("Bphen_Ag",       "Bphen + Ag (chelate)",         E.get("Bphen_Ag")),
    ("B3PyMPM_Ag",     "B3PyMPM + Ag (pyridyl N)",     E.get("B3PyMPM_Ag")),
    ("pyridine_Ag",    "Pyridine + Ag (site ref)",     E.get("pyridine_Ag")),
    ("Mo3O9_Ag",       "Mo3O9 + Ag (MoO3 model)",      E.get("Mo3O9_Ag")),
    ("Mo3O8_Ag",       "Mo3O8 + Ag (O-vacancy MoOx)",  E.get("Mo3O8_Ag")),
    ("LiF32_Ag",       "LiF(001) cluster + Ag",        E.get("LiF32_Ag")),
    ("Cs2CO3_Ag",      "Cs2CO3 + Ag",                  E.get("Cs2CO3_Ag")),
    ("F4TCNQ_Ag",      "F4TCNQ + Ag (nitrile N)",      E.get("F4TCNQ_Ag")),
    ("Liq_Ag",         "Liq + Ag",                     E.get("Liq_Ag")),
    ("TPA_Ag",         "Triphenylamine + Ag (TAPC proxy)", E.get("TPA_Ag")),
    ("PhCz_Ag",        "N-Ph-carbazole + Ag (TCTA proxy)", E.get("PhCz_Ag")),
    ("benzene_Ag",     "Benzene + Ag (pi ref)",        E.get("benzene_Ag")),
    ("HATCN_Zn",       "HATCN + Zn atom",              E.get("Zn_on_HATCN")),
    ("PhCN_ZnS",       "PhCN + ZnS molecule",          E.get("PhCN_ZnS")),
    ("Ag2",            "Ag2 dimer (reference)",        E.get("Ag2")),
]

entries = []
for tag, name, eb in SYSTEMS:
    p = os.path.join(RUNS, tag, "xtbopt.xyz")
    if not os.path.exists(p):
        p = os.path.join(RUNS, tag, "in.xyz")
    if not os.path.exists(p):
        print("skip", tag); continue
    xyz = open(p).read()
    entries.append({"tag": tag, "name": name, "eb": eb, "xyz": xyz})

html = """<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Ag Seed Layer DFT Screening — 3D Structures</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
<style>
body{font-family:'Segoe UI',sans-serif;margin:0;display:flex;height:100vh;background:#1e1e28;color:#eee}
#side{width:320px;overflow-y:auto;padding:12px;background:#26263a;box-sizing:border-box}
#view{flex:1;position:relative}
h2{font-size:15px;margin:4px 0 10px}
.btn{display:block;width:100%;text-align:left;margin:3px 0;padding:8px 10px;border:0;border-radius:8px;
 background:#33334d;color:#ddd;cursor:pointer;font-size:12.5px}
.btn:hover{background:#44446a}.btn.on{background:#5a5aa0;color:#fff}
.eb{float:right;font-weight:600;color:#9fd49f}.eb.weak{color:#e0a0a0}.eb.mid{color:#e6d38f}
#info{position:absolute;top:10px;left:10px;background:rgba(20,20,30,.85);padding:10px 14px;border-radius:10px;
 font-size:13px;z-index:5;max-width:420px}
small{color:#aaa}
</style></head><body>
<div id="side">
<h2>Ag seed-layer DFT screening<br><small>PBE-D3(BJ)/def2-SVP, CP-corrected E<sub>b</sub>(Ag)</small></h2>
__BUTTONS__
<p style="font-size:11.5px;color:#999">기준: Ag–Ag(Ag₂) 1.86 eV · Ag 벌크 응집 2.95 eV<br>
초록 ≥0.8 · 노랑 0.4–0.8 · 빨강 &lt;0.4 eV<br>클릭-드래그 회전, 휠 확대. 은색 큰 원자 = Ag</p>
</div>
<div id="view"><div id="info"></div></div>
<script>
const DATA = __DATA__;
let viewer = $3Dmol.createViewer("view",{backgroundColor:"#1e1e28"});
function show(i){
  document.querySelectorAll('.btn').forEach((b,j)=>b.classList.toggle('on',i===j));
  const d = DATA[i];
  viewer.clear();
  viewer.addModel(d.xyz,"xyz");
  viewer.setStyle({},{stick:{radius:0.12},sphere:{scale:0.25}});
  viewer.setStyle({elem:"Ag"},{sphere:{scale:0.55,color:"#c0c0c0"}});
  viewer.setStyle({elem:"Zn"},{sphere:{scale:0.5,color:"#7aa8c0"},stick:{radius:0.12}});
  viewer.setStyle({elem:"Cs"},{sphere:{scale:0.6,color:"#8f5fbf"}});
  viewer.zoomTo(); viewer.render();
  const ebtxt = d.eb==null ? "" : `E<sub>b</sub>(Ag) = <b>${d.eb.toFixed(2)} eV</b>`;
  document.getElementById("info").innerHTML = `<b>${d.name}</b><br>${ebtxt}`;
}
show(0);
</script></body></html>"""

btns = []
for i, e in enumerate(entries):
    if e["eb"] is None:
        ebs = ""
    else:
        cls = "eb" if e["eb"] >= 0.8 else ("eb mid" if e["eb"] >= 0.4 else "eb weak")
        ebs = f'<span class="{cls}">{e["eb"]:.2f} eV</span>'
    btns.append(f'<button class="btn" onclick="show({i})">{e["name"]} {ebs}</button>')

html = html.replace("__BUTTONS__", "\n".join(btns)).replace("__DATA__", json.dumps(
    [{"tag": e["tag"], "name": e["name"], "eb": e["eb"], "xyz": e["xyz"]} for e in entries]))

out = os.path.join(BASE, "structures_viewer.html")
open(out, "w", encoding="utf-8").write(html)
print("wrote", out, "with", len(entries), "structures")
