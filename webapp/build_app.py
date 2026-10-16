"""Build webapp/index.html from data/populations.yaml, energy_curves.json, fem_shapes.json and core.js.

    python webapp/build_app.py
"""
import json, pathlib, yaml
H = pathlib.Path(__file__).resolve().parent; R = H.parent
pop = yaml.safe_load((R / "data" / "populations.yaml").read_text(encoding="utf-8"))
NAMES = {"neutrophil": "Neutrophil", "lymphocyte": "Lymphocyte", "monocyte": "Monocyte", "eosinophil": "Eosinophil",
         "basophil": "Basophil", "patient_breast_CTC": "Breast CTC (patient-derived)"}
def m(s): return s["value"] if s["dist"] == "fixed" else s["mean"]
phen = [dict(id=k, name=NAMES.get(k, k), group=p["group"], D=round(m(p["D"]), 2), E_app=round(m(p["E_app"]), 3),
             gamma_pN_um=round(1000 * m(p["gamma"]), 1), tau=round(m(p["tau"]), 4), an_a=round(m(p["an_a"]), 3),
             source=f"E_app: {p['E_app']['status']}; diameter: {p['D']['status']}; γ: {p['gamma']['status']}.")
        for k, p in pop["populations"].items()]
curves = json.loads((H / "energy_curves.json").read_text())
DEV = {"membrane_PC": ("Track-etched polycarbonate membrane (10 µm thick)", "membrane"),
       "channel_PDMS": ("PDMS channel constriction", "channel")}
devices = {k: dict(name=DEV[k][0], type=DEV[k][1], profile=v["profile"],
                   u=[round(x, 5) for x in v["u"]], g=[round(x, 7) for x in v["g"]]) for k, v in curves.items()}
data = dict(K_E=pop["model"]["K_E"], phenotypes=phen, devices=devices)
core = (H / "core.js").read_text().replace("if (typeof module !== 'undefined') module.exports", "// node export:")
html = ((H / "template.html").read_text(encoding="utf-8").replace("__DATA__", json.dumps(data, ensure_ascii=False))
        .replace("__SHAPES__", (H / "fem_shapes.json").read_text()).replace("__CORE__", core))
(H / "index.html").write_text(html, encoding="utf-8")
print(f"index.html: {len(html) // 1024} kB; {len(phen)} phenotypes; {len(devices)} devices")
