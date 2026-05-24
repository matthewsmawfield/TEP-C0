"""BBN nuclear-network resolver for TEP-C0.

This module uses the installed ``BBN`` package, which integrates a stiff
Peebles/Kolb-Turner reaction network.  It is not a PArthENoPE replacement, but
it is a real network calculation and is therefore the correct resolver layer
between the C0 pipeline and future PArthENoPE/AlterBBN-grade engines.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import importlib
from typing import Any, Union

import numpy as np


C_KM_S = 299_792.458
T_CMB0_MEV = 2.7255 * 8.617333262145e-11


@dataclass(frozen=True)
class BBNInputs:
    eta: float = 6.1e-10
    h0_km_s_mpc: float = 70.0
    sigma_0: float = 0.0
    epsilon_t: float = 0.0
    n_nu: float = 3.0
    tau_n_s: float = 878.4
    n_step: int = 256
    rtol: float = 1e-4
    atol: float = 1e-10


def tep_gamma_from_temperature(T_mev: Union[float, np.ndarray], inputs: BBNInputs) -> Union[float, np.ndarray]:
    """Return the TEP expansion multiplier at photon temperature ``T``."""
    if inputs.sigma_0 == 0.0 or inputs.epsilon_t == 0.0:
        return np.ones_like(np.asarray(T_mev, dtype=float))
    one_plus_z = np.maximum(np.asarray(T_mev, dtype=float) / T_CMB0_MEV, 1.0)
    exponent = inputs.epsilon_t * inputs.sigma_0 * C_KM_S / inputs.h0_km_s_mpc
    return np.exp(exponent * np.log(one_plus_z))


@contextmanager
def _patched_bbn_expansion(inputs: BBNInputs):
    """Temporarily patch the BBN package expansion history with TEP H(T)."""
    try:
        if not hasattr(np, "int"):
            # The published package still uses np.int at import time.
            np.int = int  # type: ignore[attr-defined]

        bbn_core = importlib.import_module("BBN.BBN")
        expansion_module = importlib.import_module("BBN.expansion")
    except ImportError as e:
        raise RuntimeError(
            "BBN package not installed. Install with: pip install BBN\n"
            "Or use the fallback analytic BBN calculator."
        ) from e

    original_expansion = bbn_core.expansion

    def tep_expansion(T0, T1, N_nu=3, n_step=256):
        expansion_module.a_nu = expansion_module.rad_const * 0.875 * N_nu
        temperature_grid = np.geomspace(T0, T1, n_step)

        def expansion_eq(T, y):
            T_nu = y[0]
            E_nu = expansion_module.a_nu * T_nu**4
            E_r = expansion_module.rad_const * T**4
            P_r = E_r / 3.0
            c_r = 4.0 * E_r
            E_e, P_e, c_e = expansion_module.electron_gas(T)
            E = E_e + E_r
            P = P_e + P_r
            heat_capacity = (c_e + c_r) / T
            H_lcdm = expansion_module.G8P3 * np.sqrt(E + E_nu)
            H_tep = H_lcdm * tep_gamma_from_temperature(T, inputs)
            dy = heat_capacity / (E + P) / 3.0
            return [dy * T_nu, -dy / H_tep]

        solution = expansion_module.solve_ivp(
            expansion_eq,
            temperature_grid[[0, -1]],
            [T0, 0.0],
            t_eval=temperature_grid,
        )
        return solution.t, solution.y[0], solution.y[1]

    bbn_core.expansion = tep_expansion
    try:
        yield bbn_core
    finally:
        bbn_core.expansion = original_expansion


def _abundance_ratios(element_names: list[str], final_mass_fractions: np.ndarray) -> dict[str, float]:
    values = {name: float(value) for name, value in zip(element_names, final_mass_fractions)}
    hydrogen_mass_fraction = max(values.get("proton", 0.0), np.finfo(float).tiny)
    deuterium = values.get("deutron", 0.0)
    helium3 = values.get("helium3", 0.0)
    helium4 = values.get("helium4", 0.0)
    lithium7_total = values.get("lithium7", 0.0) + values.get("beryllium7", 0.0)
    return {
        "Y_p": helium4,
        "D_H": (deuterium / 2.0) / hydrogen_mass_fraction,
        "He3_H": (helium3 / 3.0) / hydrogen_mass_fraction,
        "Li7_H": (lithium7_total / 7.0) / hydrogen_mass_fraction,
        "hydrogen_mass_fraction": hydrogen_mass_fraction,
    }


def run_network(inputs: BBNInputs) -> dict[str, Any]:
    """Run the BBN nuclear network and return final light-element abundances."""
    with _patched_bbn_expansion(inputs) as bbn_core:
        bbn_core.initialize(N_nu=inputs.n_nu, tau_n=inputs.tau_n_s)
        element_names = list(bbn_core.e_name)
        temperature, mass_fractions = bbn_core.BBN(
            inputs.eta,
            element_names,
            n_step=inputs.n_step,
            rtol=inputs.rtol,
            atol=inputs.atol,
        )
    final_mass_fractions = mass_fractions[:, -1]
    return {
        "inputs": inputs.__dict__,
        "temperature_final_mev": float(temperature[-1]),
        "element_mass_fractions": {
            name: float(value) for name, value in zip(element_names, final_mass_fractions)
        },
        "abundances": _abundance_ratios(element_names, final_mass_fractions),
        "network": {
            "engine": "python-BBN",
            "equations": "Peebles/Kolb-Turner stiff nuclear network",
            "tep_expansion_multiplier": "H_TEP(T)=H_LCDM(T)*Gamma_TEP(T)",
        },
    }


def observational_covariance(registry_rows: list[dict[str, Any]]) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Build a diagonal observational covariance vector from the ingested registry."""
    mapping = {
        "helium_4_mass_fraction_Yp": "Y_p",
        "deuterium_to_hydrogen_DH": "D_H",
        "helium_3_to_hydrogen_He3H": "He3_H",
        "lithium_7_to_hydrogen_Li7H": "Li7_H",
    }
    labels: list[str] = []
    values: list[float] = []
    variances: list[float] = []
    for row in registry_rows:
        label = mapping.get(str(row.get("quantity")))
        if not label:
            continue
        obs_sigma = float(row["observed_sigma"])
        pred_sigma = float(row.get("prediction_sigma", 0.0))
        labels.append(label)
        values.append(float(row["observed_value"]))
        variances.append(max(obs_sigma**2 + pred_sigma**2, np.finfo(float).tiny))
    return labels, np.asarray(values, dtype=float), np.diag(variances)


def chi2_against_registry(abundances: dict[str, float], registry_rows: list[dict[str, Any]]) -> dict[str, Any]:
    labels, observed, covariance = observational_covariance(registry_rows)
    predicted = np.asarray([abundances[label] for label in labels], dtype=float)
    residual = predicted - observed
    variance = np.maximum(np.diag(covariance), np.finfo(float).tiny)
    chi2_terms = residual**2 / variance
    chi2 = float(np.sum(chi2_terms))
    return {
        "labels": labels,
        "observed": observed.tolist(),
        "predicted": predicted.tolist(),
        "covariance": covariance.tolist(),
        "residual": residual.tolist(),
        "chi2_terms": {label: float(value) for label, value in zip(labels, chi2_terms)},
        "chi2": chi2,
        "dof": len(labels),
        "reduced_chi2": float(chi2 / max(len(labels), 1)),
    }
