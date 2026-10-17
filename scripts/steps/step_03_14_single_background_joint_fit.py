#!/usr/bin/env python3
"""
Step 03.14: Single-background joint fit (Pantheon+ x DESI BAO x CMB theta*)
==========================================================================
Purpose
-------
Step 03.05 couples Pantheon+ and the CMB, but evaluates them on two different
backgrounds: the SN term uses the zero-Lambda transport kernel (Omega_m = 1,
epsilon_T ~ 0.88), while the CMB term uses hi_class with omega_cdm = 0.12 and a
drift suppressed to ~1e-6. Step 05.11 records the resulting 5.3 sigma branch
tension as an open item. This step closes that gap at the level of geometry:
every probe is evaluated on ONE expansion history per model.

Probes (all on the same background)
-----------------------------------
  SN  : Pantheon+SH0ES, z_HD > 0.01, non-calibrators, full stat+sys covariance,
        absolute magnitude M profiled analytically.
  BAO : DESI DR2 Gaussian compression (13 points, full covariance); DESI DR1
        (2024) as a cross-check. Official files: data/raw/desi_bao/.
  CMB : Planck 2018 100 theta* = 1.04110 +/- 0.00031 at z* = 1089.92.

Ruler treatment (deliberately generous to all models)
-----------------------------------------------------
The early-universe rulers are not refitted. Only their ratio is fixed,
r_s*/r_d = 144.43/147.09 (Planck 2018), which TEP preserves by conformal
closure (TEP-HC: r_s^TEP/r_s^LCDM = 0.999994). The absolute scale enters through
one free nuisance alpha = c/(H0 r_d), shared by BAO and CMB. With M profiled in
the SN term, H0 and r_d never enter separately, so the test is purely the SHAPE
of the distance-redshift relation from z = 0.01 to z = 1090.

Models
------
  LCDM   : flat, E^2 = Om(1+z)^3 + Or(1+z)^4 + (1-Om-Or)           [Om, alpha]
  M1     : corpus zero-Lambda kernel (core.cosmology.TEPCosmology with
           Omega_m = 1) plus radiation:
           dD_C/dz = (c/H0) gamma(z) / E_EdS(z),
           gamma = 1 + eps ln(1+z) exp(-(z/z_T)^n)        [eps, z_T, alpha] (n=2)
  M1n    : as M1 with n_T free                             [eps, z_T, n, alpha]
  Native : Paper 18 / hi_class branch, H = M(z) H_LCDM with
           core.cosmology.jordan_frame_M                   [Om, eps, z_T, alpha]

Output: results/step_03_14_single_background_joint_fit.json (+ .md)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve, cholesky, solve_triangular
from scipy.optimize import differential_evolution, minimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from core.cosmology import jordan_frame_M  # noqa: E402

C_KMS = 299792.458
Z_STAR = 1089.92
THETA100_OBS, THETA100_ERR = 1.04110, 0.00031
RS_STAR_OVER_RD = 144.43 / 147.09
RD_REF = 147.09
H_REF = 0.674
OMEGA_R_H2 = 4.18e-5          # photons + 3.046 massless nu
OR_REF = OMEGA_R_H2 / H_REF**2

# Corpus reference values
STEP0305 = dict(eps=0.8826330465353405, zT=4.2413196683109975, n=2.0)
STEP0305_SN_EFF = 1.0  # s_host/s_host_ref at posterior rho_half ~ ref; reported below

# ----------------------------------------------------------------- data
def load_pantheon():
    raw = ROOT / "data" / "raw"
    with open(raw / "pantheon_plus_shoes.dat") as fh:
        hdr = fh.readline().split()
        rows = [ln.split() for ln in fh if ln.strip()]
    col = {k: i for i, k in enumerate(hdr)}
    zhd = np.array([float(r[col["zHD"]]) for r in rows])
    zhel = np.array([float(r[col["zHEL"]]) for r in rows])
    mb = np.array([float(r[col["m_b_corr"]]) for r in rows])
    cal = np.array([int(r[col["IS_CALIBRATOR"]]) for r in rows])
    with open(raw / "Pantheon+SH0ES.cov") as fh:
        n = int(fh.readline())
        cov = np.array([float(x) for x in fh.read().split()]).reshape(n, n)
    assert n == len(zhd), (n, len(zhd))
    m = (zhd > 0.01) & (cal == 0)
    cov = cov[np.ix_(m, m)]
    L = cholesky(cov, lower=True)
    one = solve_triangular(L, np.ones(m.sum()), lower=True)
    return dict(z=zhd[m], zhel=zhel[m], mb=mb[m], L=L, one=one, w=float(one @ one), n=int(m.sum()))


def load_desi(tag):
    d = ROOT / "data" / "raw" / "desi_bao"
    stem = "desi_gaussian_bao_ALL_GCcomb" if tag == "DR2" else "desi_2024_gaussian_bao_ALL_GCcomb"
    z, val, kind = [], [], []
    for ln in open(d / f"{stem}_mean.txt"):
        if ln.startswith("#") or not ln.strip():
            continue
        a, b, c = ln.split()
        z.append(float(a)); val.append(float(b)); kind.append(c)
    cov = np.loadtxt(d / f"{stem}_cov.txt")
    return dict(z=np.array(z), val=np.array(val), kind=kind, cf=cho_factor(cov), n=len(z))


# ------------------------------------------------------------ backgrounds
X_GRID = np.linspace(0.0, np.log1p(Z_STAR), 6001)
Z_GRID = np.expm1(X_GRID)


def _gamma_kernel(z, eps, zT, n):
    return np.maximum(1.0 + eps * np.log1p(z) * np.exp(-(z / zT) ** n), 0.1)


def background(model, p):
    """Return (z_grid, d_C, e_eff) with d_C = D_C H0/c and e_eff = H/H0."""
    z = Z_GRID
    zp = 1.0 + z
    if model == "LCDM":
        om = p[0]
        e = np.sqrt(om * zp**3 + OR_REF * zp**4 + (1 - om - OR_REF))
        g = np.ones_like(z)
    elif model in ("M1", "M1n"):
        eps, zT = p[0], p[1]
        n = p[2] if model == "M1n" else 2.0
        e = np.sqrt((1 - OR_REF) * zp**3 + OR_REF * zp**4)
        g = _gamma_kernel(z, eps, zT, n)
    elif model == "Native":
        om, eps, zT = p[0], p[1], p[2]
        e_l = np.sqrt(om * zp**3 + OR_REF * zp**4 + (1 - om - OR_REF))
        M = np.asarray(jordan_frame_M(np.maximum(z, 1e-8), eps, zT, 2.0))
        e = e_l * M
        g = np.ones_like(z)
    else:
        raise ValueError(model)
    integrand = g / e * zp           # dz = (1+z) dx
    dC = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(X_GRID))])
    e_eff = e / g                    # c/H_eff = (c/H0) g/e  (line-of-sight ruler D_H)
    return z, dC, e_eff


def interp(zq, z, f):
    return np.interp(np.log1p(zq), np.log1p(z), f)


# --------------------------------------------------------------- chi2
def chi2_sn(sn, z, dC):
    dM = interp(sn["z"], z, dC)
    mu_shape = 5 * np.log10((1 + sn["zhel"]) * dM)   # + const absorbed in M
    y = solve_triangular(sn["L"], sn["mb"] - mu_shape, lower=True)
    off = float(sn["one"] @ y) / sn["w"]
    r = y - off * sn["one"]
    return float(r @ r)


def bao_theory(bao, z, dC, e_eff, alpha):
    dM = alpha * interp(bao["z"], z, dC)
    dH = alpha / interp(bao["z"], z, e_eff)
    out = np.empty(bao["n"])
    for i, k in enumerate(bao["kind"]):
        if k == "DM_over_rs":
            out[i] = dM[i]
        elif k == "DH_over_rs":
            out[i] = dH[i]
        elif k == "DV_over_rs":
            out[i] = (bao["z"][i] * dM[i] ** 2 * dH[i]) ** (1 / 3)
        else:
            raise ValueError(k)
    return out


def chi2_bao(bao, z, dC, e_eff, alpha):
    r = bao["val"] - bao_theory(bao, z, dC, e_eff, alpha)
    return float(r @ cho_solve(bao["cf"], r))


def theta100(dC, alpha):
    return 100 * RS_STAR_OVER_RD / (alpha * dC[-1])


def chi2_cmb(dC, alpha):
    return ((theta100(dC, alpha) - THETA100_OBS) / THETA100_ERR) ** 2


# ------------------------------------------------------------- fitting
SPECS = {
    "LCDM":   dict(names=["Om"], bounds=[(0.05, 0.95)]),
    "M1":     dict(names=["eps", "zT"], bounds=[(0.0, 3.0), (0.3, 15.0)]),
    "M1n":    dict(names=["eps", "zT", "n"], bounds=[(0.0, 3.0), (0.3, 15.0), (0.5, 6.0)]),
    "Native": dict(names=["Om", "eps", "zT"], bounds=[(0.05, 0.95), (-0.3, 0.3), (0.3, 15.0)]),
}
ALPHA_B = (15.0, 45.0)


def total(model, p, sn, bao, use, alpha):
    try:
        z, dC, e = background(model, p)
    except (ValueError, FloatingPointError):
        return 1e12, {}
    parts = {}
    if "SN" in use:
        parts["SN"] = chi2_sn(sn, z, dC)
    if "BAO" in use:
        parts["BAO"] = chi2_bao(bao, z, dC, e, alpha)
    if "CMB" in use:
        parts["CMB"] = chi2_cmb(dC, alpha)
    tot = sum(parts.values())
    return (tot if np.isfinite(tot) else 1e12), parts


def fit(model, sn, bao, use, seed=1):
    spec = SPECS[model]
    need_alpha = ("BAO" in use) or ("CMB" in use)
    bounds = spec["bounds"] + ([ALPHA_B] if need_alpha else [])

    def f(v):
        p = v[: len(spec["names"])]
        a = v[-1] if need_alpha else 30.0
        return total(model, p, sn, bao, use, a)[0]

    de = differential_evolution(f, bounds, seed=seed, tol=1e-7, maxiter=300, popsize=20, polish=False)
    res = minimize(f, de.x, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-6, maxiter=20000))
    v = res.x if res.fun <= de.fun else de.x
    p = v[: len(spec["names"])]
    a = float(v[-1]) if need_alpha else None
    chi2, parts = total(model, p, sn, bao, use, a if a else 30.0)
    out = dict(params=dict(zip(spec["names"], map(float, p))), chi2=float(chi2),
               parts={k: float(x) for k, x in parts.items()}, k=len(bounds))
    if a is not None:
        z, dC, e = background(model, p)
        out["alpha_c_over_H0rd"] = a
        out["H0_if_rd_147.09"] = C_KMS / (a * RD_REF)
        out["theta100"] = float(theta100(dC, a))
        out["theta_pull_sigma"] = float((theta100(dC, a) - THETA100_OBS) / THETA100_ERR)
        out["DM_zstar_over_LCDMref"] = float(dC[-1] / background("LCDM", [0.315])[1][-1])
    return out


def probe_at(model, p, sn, bao):
    """chi2 per probe at fixed shape params, alpha optimised on BAO+CMB."""
    z, dC, e = background(model, p)
    best = minimize(lambda a: chi2_bao(bao, z, dC, e, a[0]) + chi2_cmb(dC, a[0]), [30.0],
                    method="Nelder-Mead").x[0]
    a_bao = minimize(lambda a: chi2_bao(bao, z, dC, e, a[0]), [30.0], method="Nelder-Mead").x[0]
    return dict(SN=chi2_sn(sn, z, dC),
                BAO_alpha_from_BAO_only=chi2_bao(bao, z, dC, e, a_bao),
                CMB_with_BAO_alpha=chi2_cmb(dC, a_bao),
                theta100_with_BAO_alpha=float(theta100(dC, a_bao)),
                BAO_plus_CMB_joint_alpha=chi2_bao(bao, z, dC, e, best) + chi2_cmb(dC, best))


def apparent_w(model, p, om_eff=0.315):
    """w(z) an observer would infer if this history were LCDM-like with Om_eff."""
    z, dC, e = background(model, p)
    zz = np.linspace(0.02, 2.5, 200)
    E = interp(zz, z, e)
    rho = E**2 - om_eff * (1 + zz) ** 3
    lnr = np.log(np.where(rho > 0, rho, np.nan))
    w = -1 + np.gradient(lnr, np.log1p(zz)) / 3
    sel = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5]
    return {f"{s}": (None if not np.isfinite(np.interp(s, zz, w)) else float(np.interp(s, zz, w))) for s in sel}, \
        (float(zz[np.argmax(~(rho > 0))]) if np.any(rho <= 0) else None)


def main():
    sn = load_pantheon()
    out = dict(step="step_03_14_single_background_joint_fit",
               data=dict(SN=f"Pantheon+SH0ES, N={sn['n']} (z_HD>0.01, non-calibrators), full cov",
                         CMB=f"Planck 2018 100theta*={THETA100_OBS}+/-{THETA100_ERR}, z*={Z_STAR}",
                         ruler=f"r_s*/r_d={RS_STAR_OVER_RD:.5f} fixed; alpha=c/(H0 r_d) free, shared"),
               fits={})
    for tag in ("DR2", "DR1"):
        bao = load_desi(tag)
        res = {}
        for use in (("SN",), ("SN", "BAO"), ("SN", "BAO", "CMB")):
            key = "+".join(use)
            res[key] = {m: fit(m, sn, bao, use) for m in SPECS}
            base = res[key]["LCDM"]
            for m in SPECS:
                r = res[key][m]
                r["dchi2_vs_LCDM"] = r["chi2"] - base["chi2"]
                r["dAIC_vs_LCDM"] = r["dchi2_vs_LCDM"] + 2 * (r["k"] - base["k"])
            print(f"[{tag}] {key}: " + ", ".join(f"{m} chi2={res[key][m]['chi2']:.2f}" for m in SPECS), flush=True)
        res["corpus_step0305_point"] = probe_at("M1", [STEP0305["eps"], STEP0305["zT"]], sn, bao)
        res["corpus_step0305_point_LCDM_ref"] = probe_at("LCDM", [res["SN+BAO+CMB"]["LCDM"]["params"]["Om"]], sn, bao)
        out["fits"][f"DESI_{tag}"] = res
    # apparent w for the joint best fits (DR2)
    j = out["fits"]["DESI_DR2"]["SN+BAO+CMB"]
    aw = {}
    for m in ("M1", "M1n", "Native"):
        p = list(j[m]["params"].values())
        om_eff = j[m]["params"].get("Om", j["LCDM"]["params"]["Om"])
        w, zneg = apparent_w(m, p, om_eff)
        aw[m] = dict(Om_eff_used=om_eff, w=w, rhoDE_negative_from_z=zneg)
    out["apparent_w_DR2_joint"] = aw
    # profile of the native branch over eps (DR2, SN+BAO+CMB), Om, z_T, alpha re-optimised
    bao = load_desi("DR2")
    L = j["LCDM"]["chi2"]
    prof = {}
    for eps in (0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04):
        best = None
        for zt0 in (0.35, 0.5, 1.0, 2.0, 5.0):
            r = minimize(lambda v: total("Native", [v[0], eps, v[1]], sn, bao, ("SN", "BAO", "CMB"), v[2])[0],
                         [0.3, zt0, 30.0], method="Nelder-Mead", options=dict(maxiter=8000, xatol=1e-7, fatol=1e-6))
            if best is None or r.fun < best.fun:
                best = r
        prof[f"{eps}"] = dict(dchi2_vs_LCDM=float(best.fun - L), Om=float(best.x[0]), zT=float(best.x[1]))
    out["native_eps_profile_DR2_joint"] = prof
    out["native_corpus_point_eps0.0043_zT5_DR2_joint_dchi2"] = float(
        minimize(lambda v: total("Native", [v[0], 0.0043, 5.0], sn, bao, ("SN", "BAO", "CMB"), v[1])[0],
                 [0.3, 30.0], method="Nelder-Mead").fun - L)
    path = ROOT / "results" / "step_03_14_single_background_joint_fit.json"
    path.write_text(json.dumps(out, indent=1))
    print("wrote", path)


if __name__ == "__main__":
    main()
