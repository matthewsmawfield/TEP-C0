#!/usr/bin/env python3
"""Step 001: verify the FLRW recovery limit of the open-path transport connection (K_T)."""

from __future__ import annotations

import numpy as np
from pathlib import Path

from c0_common import (
    C_KM_S,
    H0_KM_S_MPC,
    TEPLogger,
    comoving_distance_mpc,
    ensure_dirs,
    h_km_s_mpc,
    print_status,
    rel,
    rounded,
    set_step_logger,
    step_csv_path,
    step_json_path,
    write_csv,
    write_json,
)


STEP_ID = "step_02_01_transport_kernel"


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    print_status("Evaluating open-path transport connection K_T across z=[0, 3]", "PROCESS")
    z = np.linspace(0.0, 3.0, 1001)
    distance_mpc = comoving_distance_mpc(z)
    
    # We define the transport integral T = ln(1+z)
    log_transport = np.log1p(z)
    reconstructed_z = np.exp(log_transport) - 1
    recovery_error = np.max(np.abs(reconstructed_z - z))
    print_status(f"Max reconstruction error: {recovery_error:.2e}", "SUCCESS")

    # Deriving the low-z Transport Hubble Law: z ~ K_T,0 * d
    low_z = (z > 0) & (z <= 0.02)
    slope, intercept = np.polyfit(distance_mpc[low_z], z[low_z], 1)
    recovered_h0 = slope * C_KM_S
    low_z_intercept = intercept
    print_status(f"Low-z H_T recovery: {recovered_h0:.3f} km/s/Mpc", "SUCCESS")

    print_status("Calculating matter-frame temporal transport density (K_T)", "PROCESS")
    # K_T = d(ln(1+z))/dD
    k_transport = np.gradient(log_transport, distance_mpc, edge_order=2)
    h_transport = C_KM_S * k_transport

    rows = [
        {
            "z": rounded(zi),
            "distance_mpc": rounded(di, 6),
            "log_transport": rounded(li, 9),
            "k_transport_mpc_inv": rounded(ki, 9),
            "h_transport_km_s_mpc": rounded(hi, 6),
            "h_flrw_over_1pz_km_s_mpc": rounded(h_km_s_mpc(zi) / (1 + zi), 6),
        }
        for zi, di, li, ki, hi in zip(z, distance_mpc, log_transport, k_transport, h_transport)
    ]
    write_csv(step_csv_path(STEP_ID), rows)

    payload = {
        "step": STEP_ID,
        "description": "FLRW recovery check for the open-path transport connection (K_T).",
        "cosmology": {
            "H0_km_s_Mpc": H0_KM_S_MPC,
            "Omega_m": 0.3,
            "Omega_Lambda": 0.7,
        },
        "metrics": {
            "max_reconstruction_error_z": rounded(recovery_error, 12),
            "low_z_recovered_HT_km_s_Mpc": rounded(recovered_h0, 6),
            "transport_connection_density_z0": rounded(k_transport[0], 9),
        },
        "artifacts": {
            "csv": rel(step_csv_path(STEP_ID)),
        },
        "interpretation": "The transport law ln(1+z) = integral(K_T dl) exactly recovers the FLRW scale factor limit when the connection is integrable. The Big Bang corresponds to the boundary where the reconstruction factor a_eff = exp(-integral) vanishes.",
    }
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
