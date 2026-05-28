from cobaya.model import get_model
import yaml
import sys
sys.path.insert(0, "core")
with open('results/outputs/tep_cobaya_sne.updated.yaml', 'r') as f:
    info = yaml.safe_load(f)
info["output"] = None
# Remove planck
del info["likelihood"]["planck_2018_highl_plik.TTTEEE"]
if "planck_2018_lowl.TT" in info["likelihood"]: del info["likelihood"]["planck_2018_lowl.TT"]
if "planck_2018_lowl.EE" in info["likelihood"]: del info["likelihood"]["planck_2018_lowl.EE"]

model = get_model(info)
point = model.prior.reference()
print("Evaluating point for SNe only")
res = model.logposterior(point)
print(f"Logpost: {res.logpost}")
print(res.loglikes)
