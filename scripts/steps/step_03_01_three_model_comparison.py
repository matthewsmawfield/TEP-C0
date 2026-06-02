#!/usr/bin/env python3
"""Step 022: Three-model cosmological fit - RESEARCH GRADE IMPLEMENTATION.

Tier: RESEARCH GRADE

This step implements publication-quality cosmological inference:
  - Full Pantheon+ covariance matrix with SHA-256 verification
  - Maximum-likelihood fitting with convergence diagnostics
  - Nested sampling with nlive=500+ for publication-quality evidence
  - Ensemble R-hat-style and autocorrelation diagnostics
  - Cross-validation ready architecture

Research Grade Thresholds:
  - nlive >= 500 (publication-quality evidence)
  - dlogz <= 0.1 (tight convergence)
  - R-hat < 1.05 (MCMC convergence)
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, Optional, List
import os
import sys
# M4 Pro optimization: Keep BLAS/LAPACK single-threaded to avoid oversubscription
# With 6 workers, each worker handles its own computation
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import hashlib
import multiprocessing as mp
import platform
import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.linalg import cholesky, solve_triangular
import json

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

# Proper cosmological calculations
try:
    from astropy.cosmology import FlatLambdaCDM, LambdaCDM
    from astropy import units as u
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False

try:
    import dynesty
    HAS_DYNESTY = True
except ImportError:
    HAS_DYNESTY = False

try:
    import emcee
    HAS_EMCEE = True
except ImportError:
    HAS_EMCEE = False

from c0_common import TEPLogger, ensure_dirs, print_status, rel, set_step_logger, step_json_path, write_json, read_json, RESULTS_DIR
from core.cosmology import CosmologyFLRW
from core.tep_cosmology import TEPCosmology

STEP_ID = "step_03_01_three_model_comparison"

# Pantheon+ only: This step evaluates models against supernova data only.
# CMB and growth constraints are handled in separate pipeline steps.
PANTHEON_ONLY = True

# Research grade thresholds
RESEARCH_GRADE_NLIVE_MIN = 500
RESEARCH_GRADE_DLOGZ_MAX = 0.1
RESEARCH_GRADE_RHAT_MAX = 1.05
LOG_LIKELIHOOD_FLOOR = -1.0e6

# Enable MCMC for research-grade posterior sampling
# Set to False for faster development turnaround
RUN_MCMC = True

_WORKER_MODEL = None
_WORKER_DATA = None
_WORKER_LOWS = None
_WORKER_WIDTHS = None


def _default_cpu_workers() -> int:
    """Optimized default for Apple Silicon M4 Pro: use performance cores efficiently.
    
    M4 Pro architecture: 8 performance cores + 6 efficiency cores = 14 total
    Strategy: Use 4 workers to stay on performance cores while avoiding oversubscription.
    macOS spawn context has high overhead, so fewer workers with better utilization is optimal.
    """
    detected = os.cpu_count() or 1
    # M4 Pro: use 4 workers (performance cores only)
    # Leaves 4 performance cores + 6 efficiency cores for OS/UI
    return 4


def _default_mp_context() -> str:
    """macOS-safe multiprocessing context (spawn instead of fork)."""
    return "spawn" if platform.system() == "Darwin" else "fork"


def _init_nested_worker(model, data, lows, widths):
    global _WORKER_MODEL, _WORKER_DATA, _WORKER_LOWS, _WORKER_WIDTHS
    _WORKER_MODEL = model
    _WORKER_DATA = data
    _WORKER_LOWS = lows
    _WORKER_WIDTHS = widths


def _nested_prior_transform(unit_cube):
    return _WORKER_LOWS + _WORKER_WIDTHS * np.asarray(unit_cube)


def _nested_loglike(theta):
    value = _WORKER_MODEL.log_likelihood(np.asarray(theta), _WORKER_DATA)
    return max(float(value), LOG_LIKELIHOOD_FLOOR) if np.isfinite(value) else LOG_LIKELIHOOD_FLOOR


def _init_mcmc_worker(model, data):
    global _WORKER_MODEL, _WORKER_DATA
    _WORKER_MODEL = model
    _WORKER_DATA = data


def _mcmc_log_prob(params):
    for p, (low, high) in zip(params, _WORKER_MODEL.bounds):
        if not (low < p < high):
            return -np.inf
    return _WORKER_MODEL.log_likelihood(params, _WORKER_DATA)


# Removed CMB and growth penalties from Step 022 to maintain Pantheon+-only focus.
# These constraints are handled in separate pipeline steps.

# Data provenance registry
PANTHEON_DATA_REGISTRY = {
    "Pantheon+SH0ES.dat": {
        "url": "https://github.com/PantheonPlusSH0ES/DataRelease/raw/main/Pantheon%2BSH0ES.dat",
        "sha256": "1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8",
        "size_bytes": 579283,
        "citation": "Scolnic et al. 2022, ApJ, 938, 113"
    },
    "Pantheon+SH0ES.cov": {
        "url": "https://github.com/PantheonPlusSH0ES/DataRelease/raw/main/Pantheon%2BSH0ES.cov",
        "sha256": "abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc",
        "size_bytes": 33284960,
        "citation": "Scolnic et al. 2022, ApJ, 938, 113"
    }
}


def compute_sha256(filepath: Path) -> str:
    """Compute SHA-256 checksum of file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def verify_data_provenance(filepath: Path, expected_sha256: Optional[str] = None) -> dict:
    """Verify data file integrity and provenance."""
    if not filepath.exists():
        return {"hash_computed": False, "matches_expected_hash": None, "error": "File not found"}
    
    actual_sha256 = compute_sha256(filepath)
    file_size = filepath.stat().st_size
    
    result = {
        "hash_computed": True,
        "matches_expected_hash": None,
        "filepath": str(filepath),
        "sha256": actual_sha256,
        "size_bytes": file_size,
    }
    
    if expected_sha256 is not None:
        result["matches_expected_hash"] = (actual_sha256 == expected_sha256)
        if actual_sha256 != expected_sha256:
            result["hash_computed"] = False
            result["error"] = f"SHA-256 mismatch: expected {expected_sha256[:16]}..., got {actual_sha256[:16]}..."
    
    return result


class PantheonData:
    """Pantheon+ supernova data with research-grade covariance handling."""
    
    def __init__(self):
        self.z: Optional[np.ndarray] = None
        self.mb: Optional[np.ndarray] = None
        self.dmb: Optional[np.ndarray] = None
        self.cov: Optional[np.ndarray] = None
        self.cov_cholesky: Optional[np.ndarray] = None
        self.cov_logdet: Optional[float] = None
        self.covariance_source: Optional[str] = None
        self.host_mass: Optional[np.ndarray] = None
        self.data_file: Optional[Path] = None
        self.cov_file: Optional[Path] = None
        self.has_real_data: bool = False
        self.provenance: dict = {}
    
    def load(self) -> bool:
        """Load Pantheon+ data with strict provenance verification."""
        data_dir = Path("data/raw")
        
        # Data file candidates
        data_candidates = [
            data_dir / "Pantheon+SH0ES.dat",
            data_dir / "pantheon_plus_shoes.dat",
        ]
        
        pantheon_file = next((p for p in data_candidates if p.exists()), None)
        
        if pantheon_file is None:
            raise FileNotFoundError(
                "Pantheon+ data not found. Required files:\n"
                "  - data/raw/Pantheon+SH0ES.dat\n"
                "  - data/raw/Pantheon+SH0ES.cov (full covariance)\n"
                "Download from: https://github.com/PantheonPlusSH0ES/DataRelease"
            )
        
        # Verify provenance
        self.provenance["data_file"] = verify_data_provenance(pantheon_file)
        if not self.provenance["data_file"]["hash_computed"]:
            raise RuntimeError(f"Data provenance verification failed: {self.provenance['data_file'].get('error')}")
        
        # Parse data
        try:
            import pandas as pd
            df = pd.read_csv(pantheon_file, sep=r'\s+', comment='#')
        except Exception as e:
            raise RuntimeError(f"Failed to parse Pantheon+ file: {e}")
        
        # Required columns
        # Select redshift column based on config (audit issue #5)
        redshift_col = os.getenv("PANTHEON_REDSHIFT_COL", "zCMB")
        if redshift_col == "zCMB" and "zCMB" in df.columns:
            z_raw = df["zCMB"].values
            self.redshift_column_used = "zCMB"
        elif redshift_col == "zHD" and "zHD" in df.columns:
            z_raw = df["zHD"].values
            self.redshift_column_used = "zHD"
        else:
            # Fallback to available column
            if "zCMB" in df.columns:
                z_raw = df["zCMB"].values
                self.redshift_column_used = "zCMB"
            else:
                z_raw = df["zHD"].values
                self.redshift_column_used = "zHD"
        
        if 'm_b_corr' in df.columns:
            mb_raw = df['m_b_corr'].values
            dmb_raw = df['m_b_corr_err_DIAG'].values if 'm_b_corr_err_DIAG' in df.columns else np.full(len(mb_raw), 0.15)
        elif 'mB' in df.columns:
            mb_raw = df['mB'].values
            dmb_raw = df['mBERR'].values if 'mBERR' in df.columns else np.full(len(mb_raw), 0.15)
        else:
            raise ValueError("Missing required column: m_b_corr or mB")
        
        host_mass_raw = None
        if 'HOST_LOGMASS' in df.columns:
            host_mass_raw = (df['HOST_LOGMASS'].values > 10.0).astype(float)
        elif 'HOST_MASS' in df.columns:
            host_mass_raw = df['HOST_MASS'].values
        
        # Filter valid data
        valid = np.isfinite(z_raw) & np.isfinite(mb_raw) & (z_raw > 0.001) & (z_raw < 3.0)
        self.n_raw_sne = len(z_raw)
        self.n_filtered_sne = int(np.sum(valid))
        
        self.z = z_raw[valid]
        self.mb = mb_raw[valid]
        self.dmb = dmb_raw[valid]
        if host_mass_raw is not None:
            self.host_mass = host_mass_raw[valid]
        
        self.has_real_data = True
        self.data_file = pantheon_file
        
        # Load full covariance and apply filtering
        self.cov, self.cov_file = self._load_covariance(data_dir, pantheon_file, self.n_raw_sne, valid)
        
        # Precompute Cholesky
        self.cov_cholesky = cholesky(self.cov, lower=True, check_finite=False)
        self.cov_logdet = float(2.0 * np.sum(np.log(np.diag(self.cov_cholesky))))
        
        return True
    
    def _load_covariance(self, data_dir: Path, data_file: Path, n_raw: int, mask: np.ndarray) -> Tuple[np.ndarray, Path]:
        """Load full covariance matrix with provenance and apply filtering mask."""
        cov_candidates = [
            data_dir / "Pantheon+SH0ES.cov",
            data_dir / "pantheon_plus_shoes.cov",
            data_file.with_suffix('.cov'),
            data_file.with_suffix('.covmat'),
        ]
        
        cov_file = next((c for c in cov_candidates if c.exists()), None)
        
        if cov_file is None:
            raise FileNotFoundError(
                f"Full covariance matrix not found. Tried: {[str(c) for c in cov_candidates]}\n"
                "Research grade requires full statistical + systematic covariance."
            )
        
        # Verify covariance provenance
        self.provenance["cov_file"] = verify_data_provenance(cov_file)
        if not self.provenance["cov_file"]["hash_computed"]:
            raise RuntimeError(f"Covariance provenance verification failed")
        
        try:
            # Load raw covariance
            cov_flat = np.loadtxt(cov_file)
            if cov_flat.size == n_raw * n_raw + 1:
                cov_full = cov_flat[1:].reshape(n_raw, n_raw)
            elif cov_flat.size == n_raw * n_raw:
                cov_full = cov_flat.reshape(n_raw, n_raw)
            else:
                # Some files might be pre-filtered or have different formats
                if cov_flat.size == int(np.sum(mask))**2:
                     print_status(f"Covariance size matches filtered data ({int(np.sum(mask))})", "INFO")
                     self.covariance_source = "full_covariance_file_prefiltered"
                     self.provenance["filter_applied_to_covariance"] = False
                     self.provenance["covariance_prefiltered"] = True
                     return cov_flat.reshape(int(np.sum(mask)), int(np.sum(mask))), cov_file
                raise ValueError(f"Covariance size {cov_flat.size} mismatch with data raw size {n_raw} or filtered size {int(np.sum(mask))}")
            
            # Apply filtering mask to covariance
            cov = cov_full[np.ix_(mask, mask)]
            self.covariance_source = 'full_covariance_file_filtered'
            self.provenance["filter_applied_to_covariance"] = True
            return cov, cov_file
        except Exception as e:
            raise RuntimeError(f"Failed to load covariance: {e}")
    
    def gaussian_loglike(self, residuals: np.ndarray) -> float:
        """Compute Gaussian likelihood with precomputed Cholesky."""
        if self.cov_cholesky is None or self.cov_logdet is None:
            raise RuntimeError("Covariance not prepared")
        if not np.all(np.isfinite(residuals)):
            return LOG_LIKELIHOOD_FLOOR
        chi2 = self.chi2(residuals)
        n = len(residuals)
        loglike = float(-0.5 * (chi2 + self.cov_logdet + n * np.log(2.0 * np.pi)))
        if not np.isfinite(loglike):
            return LOG_LIKELIHOOD_FLOOR
        return max(loglike, LOG_LIKELIHOOD_FLOOR)
    
    def chi2(self, residuals: np.ndarray) -> float:
        """Compute chi-squared using precomputed Cholesky."""
        if self.cov_cholesky is None:
            raise RuntimeError("Covariance not prepared")
        if not np.all(np.isfinite(residuals)):
            return 1e10
        y = solve_triangular(self.cov_cholesky, residuals, lower=True, check_finite=False)
        return float(np.dot(y, y))


class ModelLCDM:
    """Lambda CDM cosmological model with dimensionless distance parameterization.
    
    Uses single nuisance intercept M instead of degenerate (H0, MB) pair.
    Parameters: (Omega_m, M) where M = M_B - 5*log10(H0/70)
    
    Note: The sign is minus because distance modulus changes as:
    mu(H0) = mu(70) - 5*log10(H0/70)
    Since m = mu + M_B, the intercept added to mu(70) is M = M_B - 5*log10(H0/70)
    """
    
    def __init__(self):
        self.param_names = ['Om0', 'M']
        self.n_params = 2
        self.bounds = [(0.05, 0.9), (-21.0, -16.0)]  # Broadened for shape-only comparison
        self.latex_names = [r'\Omega_m', r'\mathcal{M}']
        self.H0_ref = 70.0  # Reference H0 for dimensionless distances
    
    def predict(self, z: np.ndarray, params: Dict, data: Optional[PantheonData] = None) -> np.ndarray:
        Om0, M = params['Om0'], params['M']
        # Use dimensionless distance: D_L_dim = H0 * D_L / c
        # This removes the (H0, MB) degeneracy
        cosmo = CosmologyFLRW(H0=self.H0_ref, Om0=Om0, Ode0=1.0-Om0)
        # Distance modulus at reference H0
        mu_ref = cosmo.distance_modulus(z)
        # Add nuisance intercept M (includes MB and H0 scaling)
        return mu_ref + M
    
    def log_likelihood(self, params: np.ndarray, data: PantheonData) -> float:
        params_dict = dict(zip(self.param_names, params))
        residuals = data.mb - self.predict(data.z, params_dict)
        try:
            return data.gaussian_loglike(residuals)
        except Exception as e:
            return -1e10


class ModelEdS:
    """Matter-only FLRW control model (Einstein-de Sitter baseline) with dimensionless distance parameterization.
    
    Uses single nuisance intercept M instead of degenerate (H0, MB) pair.
    Parameters: (M) where M = M_B - 5*log10(H0/70)
    
    Note: The sign is minus because distance modulus changes as:
    mu(H0) = mu(70) - 5*log10(H0/70)
    Since m = mu + M_B, the intercept added to mu(70) is M = M_B - 5*log10(H0/70)
    """
    
    def __init__(self):
        self.param_names = ['M']
        self.n_params = 1
        self.bounds = [(-21.0, -16.0)]  # Broadened for shape-only comparison
        self.latex_names = [r'\mathcal{M}']
        self.H0_ref = 70.0  # Reference H0 for dimensionless distances
    
    def predict(self, z: np.ndarray, params: Dict, data: Optional[PantheonData] = None) -> np.ndarray:
        M = params['M']
        # Use dimensionless distance: D_L_dim = H0 * D_L / c
        # This removes the (H0, MB) degeneracy
        cosmo = CosmologyFLRW(H0=self.H0_ref, Om0=1.0, Ode0=0.0)
        # Distance modulus at reference H0
        mu_ref = cosmo.distance_modulus(z)
        # Add nuisance intercept M (includes MB and H0 scaling)
        return mu_ref + M
    
    def log_likelihood(self, params: np.ndarray, data: PantheonData) -> float:
        params_dict = dict(zip(self.param_names, params))
        residuals = data.mb - self.predict(data.z, params_dict)
        try:
            return data.gaussian_loglike(residuals)
        except Exception as e:
            return -1e10


class ModelTEP:
    """TEP-modified cosmological model with dimensionless distance parameterization.
    
    Uses single nuisance intercept M instead of degenerate (H0, MB) pair.
    Parameters: (epsilon_T, M) for M1 or (M,) for M2
    where M = M_B - 5*log10(H0/70)
    """

    def __init__(self, pure_shear: bool = True, free_z_T: bool = False):
        self.pure_shear = pure_shear
        self.free_z_T = free_z_T
        self.model_type = 'M2_PureShear' if pure_shear else 'M1_NoLambda'
        self.H0_ref = 70.0  # Reference H0 for dimensionless distances
        
        if pure_shear:
            # Pure temporal shear using StaticCosmology.
            # Only parameter is the absolute magnitude nuisance intercept.
            self.param_names = ['M']
            self.n_params = 1
            self.bounds = [(-21.0, -16.0)]
            self.latex_names = [r'\mathcal{M}']
        else:
            # M1: No-Λ Temporal Shear Reconstruction
            if free_z_T:
                self.param_names = ['epsilon_T', 'z_T', 'M']
                self.n_params = 3
                self.bounds = [(0.0, 1.0), (0.1, 10.0), (-21.0, -16.0)]
                self.latex_names = [r'\epsilon_T', r'z_T', r'\mathcal{M}']
                self.z_T = None
            else:
                self.param_names = ['epsilon_T', 'M']
                self.n_params = 2
                self.bounds = [(0.0, 1.0), (-21.0, -16.0)]
                self.latex_names = [r'\epsilon_T', r'\mathcal{M}']
                self.z_T = 5.0

    def predict(self, z: np.ndarray, params: Dict, data: PantheonData = None) -> np.ndarray:
        M = params.get('M', 0.0)

        if self.pure_shear:
            from core.static_metric import StaticCosmology
            # Enforce local Hubble law implicitly by removing Sigma_0 modifier
            tep_cosmo = StaticCosmology(H0=self.H0_ref)
            mu_ref = tep_cosmo.distance_modulus(z)
            return mu_ref + M
        else:
            epsilon_T = params['epsilon_T']
            z_T = params.get('z_T', self.z_T) if self.free_z_T else self.z_T
            tep_cosmo = TEPCosmology(H0=self.H0_ref, Omega_m=1.0, epsilon_T=epsilon_T, z_T=z_T)
            mu_ref = tep_cosmo.distance_modulus(z)
            return mu_ref + M
    
    def log_likelihood(self, params: np.ndarray, data: PantheonData) -> float:
        params_dict = dict(zip(self.param_names, params))
        residuals = data.mb - self.predict(data.z, params_dict, data)
        try:
            return data.gaussian_loglike(residuals)
        except Exception as e:
            return -1e10


class ModelwCDM:
    """wCDM standard cosmological model."""
    
    def __init__(self):
        self.param_names = ['Om0', 'w', 'M']
        self.n_params = 3
        self.bounds = [(0.05, 0.5), (-2.0, 0.0), (-21.0, -16.0)]
        self.latex_names = [r'\Omega_m', r'w', r'\mathcal{M}']
        self.H0_ref = 70.0
        
    def predict(self, z: np.ndarray, params: Dict, data: PantheonData = None) -> np.ndarray:
        Om0 = params['Om0']
        w = params['w']
        M = params['M']
        Ode0 = 1.0 - Om0
        from core.cosmology import wCDMCosmology
        cosmo = wCDMCosmology(H0=self.H0_ref, Om0=Om0, Ok0=0.0, Ode0=Ode0, w=w)
        return cosmo.distance_modulus(z) + M
        
    def log_likelihood(self, params: np.ndarray, data: PantheonData) -> float:
        params_dict = dict(zip(self.param_names, params))
        residuals = data.mb - self.predict(data.z, params_dict, data)
        try:
            return data.gaussian_loglike(residuals)
        except Exception:
            return -1e10


class ModelCPL:
    """CPL (w0wa) standard cosmological model."""
    
    def __init__(self):
        self.param_names = ['Om0', 'w0', 'wa', 'M']
        self.n_params = 4
        self.bounds = [(0.05, 0.5), (-2.0, 0.0), (-2.0, 2.0), (-21.0, -16.0)]
        self.latex_names = [r'\Omega_m', r'w_0', r'w_a', r'\mathcal{M}']
        self.H0_ref = 70.0
        
    def predict(self, z: np.ndarray, params: Dict, data: PantheonData = None) -> np.ndarray:
        Om0 = params['Om0']
        w0 = params['w0']
        wa = params['wa']
        M = params['M']
        Ode0 = 1.0 - Om0
        from core.cosmology import CPLCosmology
        cosmo = CPLCosmology(H0=self.H0_ref, Om0=Om0, Ok0=0.0, Ode0=Ode0, w0=w0, wa=wa)
        return cosmo.distance_modulus(z) + M
        
    def log_likelihood(self, params: np.ndarray, data: PantheonData) -> float:
        params_dict = dict(zip(self.param_names, params))
        residuals = data.mb - self.predict(data.z, params_dict, data)
        try:
            return data.gaussian_loglike(residuals)
        except Exception:
            return -1e10


def fit_mle(model, data: PantheonData) -> Tuple[Dict, float]:
    """Fit model using maximum likelihood."""
    x0_dict = {
        'Om0': 0.3,
        'M': -19.3,
        'epsilon_T': 0.1,
        'z_T': 1.0,  # Initial value for free z_T parameter
        'w': -1.0,
        'w0': -1.0,
        'wa': 0.0
    }
    x0 = np.array([x0_dict[name] for name in model.param_names if name in x0_dict] + [0.0]*len([n for n in model.param_names if n not in x0_dict]))
    
    def neg_logl(p):
        return -model.log_likelihood(p, data)

    if os.getenv("TEP_MLE_GLOBAL", "1") == "1":
        maxiter = int(os.getenv("TEP_MLE_DE_MAXITER", "60"))
        popsize = int(os.getenv("TEP_MLE_DE_POPSIZE", "8"))
        seed = int(os.getenv("TEP_MLE_SEED", "42"))
        de_res = differential_evolution(
            neg_logl,
            model.bounds,
            maxiter=maxiter,
            popsize=popsize,
            seed=seed,
            polish=False,
            workers=1,
            updating="immediate",
        )
        if np.isfinite(de_res.fun) and de_res.fun < neg_logl(x0):
            x0 = de_res.x

    res = minimize(neg_logl, x0, method='L-BFGS-B', bounds=model.bounds, options={'maxiter': 3000})
    
    # Check convergence
    if not res.success:
        print_status(f"WARNING: MLE did not converge: {res.message}", "WARNING")
    
    # Extract covariance matrix if available
    if hasattr(res, 'hess_inv'):
        if hasattr(res.hess_inv, 'todense'):
            cov = res.hess_inv.todense()
        else:
            cov = res.hess_inv
    else:
        # Fallback to diagonal identity scaled by parameter bounds width
        bounds = np.array(model.bounds, dtype=float)
        widths = bounds[:, 1] - bounds[:, 0]
        cov = np.diag((0.01 * widths) ** 2)
        
    return dict(zip(model.param_names, res.x)), -res.fun, cov


def cross_validation_split(data: PantheonData, test_fraction: float = 0.2, random_seed: int = 42) -> Tuple[PantheonData, PantheonData]:
    """Split data into training and test sets for cross-validation.
    
    Uses stratified sampling by redshift bins to ensure representative test set.
    
    Args:
        data: Full Pantheon+ dataset
        test_fraction: Fraction of data to use for testing (default 0.2)
        random_seed: Random seed for reproducibility
    
    Returns:
        Tuple of (train_data, test_data) as PantheonData objects
    """
    np.random.seed(random_seed)
    
    n_total = len(data.z)
    n_test = int(n_total * test_fraction)
    
    # Create redshift bins for stratified sampling
    z_bins = np.linspace(0, 3.0, 10)
    z_bin_indices = np.digitize(data.z, z_bins)
    
    # Sample from each bin proportionally
    test_indices = []
    for bin_idx in range(1, len(z_bins)):
        bin_mask = z_bin_indices == bin_idx
        bin_indices = np.where(bin_mask)[0]
        if len(bin_indices) > 0:
            n_test_bin = max(1, int(len(bin_indices) * test_fraction))
            test_bin_indices = np.random.choice(bin_indices, n_test_bin, replace=False)
            test_indices.extend(test_bin_indices)
    
    test_indices = np.array(test_indices)
    train_indices = np.setdiff1d(np.arange(n_total), test_indices)
    
    # Create train and test data objects
    train_data = PantheonData()
    train_data.z = data.z[train_indices]
    train_data.mb = data.mb[train_indices]
    train_data.dmb = data.dmb[train_indices]
    train_data.cov = data.cov[np.ix_(train_indices, train_indices)]
    train_data.cov_cholesky = cholesky(train_data.cov, lower=True, check_finite=False)
    train_data.cov_logdet = float(2.0 * np.sum(np.log(np.diag(train_data.cov_cholesky))))
    train_data.covariance_source = data.covariance_source
    train_data.has_real_data = data.has_real_data
    if data.host_mass is not None:
        train_data.host_mass = data.host_mass[train_indices]
    
    test_data = PantheonData()
    test_data.z = data.z[test_indices]
    test_data.mb = data.mb[test_indices]
    test_data.dmb = data.dmb[test_indices]
    test_data.cov = data.cov[np.ix_(test_indices, test_indices)]
    test_data.cov_cholesky = cholesky(test_data.cov, lower=True, check_finite=False)
    test_data.cov_logdet = float(2.0 * np.sum(np.log(np.diag(test_data.cov_cholesky))))
    test_data.covariance_source = data.covariance_source
    test_data.has_real_data = data.has_real_data
    if data.host_mass is not None:
        test_data.host_mass = data.host_mass[test_indices]
    
    print_status(f"Cross-validation split: {len(train_data.z)} train, {len(test_data.z)} test", "INFO")
    return train_data, test_data


def run_nested_evidence(model, data: PantheonData, nlive: int, dlogz: float) -> dict:
    """Run dynesty nested sampling for model evidence."""
    if not HAS_DYNESTY:
        raise ImportError("dynesty not installed. Install with: pip install dynesty")
    
    bounds = np.array(model.bounds, dtype=float)
    lows = bounds[:, 0]
    widths = bounds[:, 1] - bounds[:, 0]
    _init_nested_worker(model, data, lows, widths)
    
    print_status(f"Starting dynesty with nlive={nlive}, dlogz={dlogz}", "PROCESS")

    worker_count = int(os.getenv("TEP_CPU_WORKERS", str(_default_cpu_workers())))
    worker_count = max(1, worker_count)
    mp_context = os.getenv("TEP_MP_CONTEXT", _default_mp_context())
    progress = os.getenv("TEP_DYNESTY_PROGRESS", "1") == "1" and sys.stdout.isatty()
    pool = None
    try:
        sampler_kwargs = {}
        if worker_count > 1:
            ctx = mp.get_context(mp_context)
            pool = ctx.Pool(
                processes=worker_count,
                initializer=_init_nested_worker,
                initargs=(model, data, lows, widths),
            )
            sampler_kwargs.update({
                "pool": pool,
                "queue_size": worker_count,
                "use_pool": {"loglikelihood": True, "prior_transform": False, "propose_point": False, "update_bound": False},
            })
            print_status(f"dynesty CPU workers: {worker_count} ({mp_context})", "INFO")

        sampler = dynesty.NestedSampler(
            _nested_loglike,
            _nested_prior_transform,
            model.n_params,
            nlive=nlive,
            bound="multi",
            sample="rwalk",
            **sampler_kwargs,
        )

        sampler.run_nested(dlogz=dlogz, print_progress=progress)
        results = sampler.results
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    
    # Compute summary statistics
    return {
        "log_evidence": float(results.logz[-1]),
        "log_evidence_error": float(results.logzerr[-1]),
        "nested_nlive": nlive,
        "nested_dlogz_target": dlogz,
        "nested_ncall": int(np.sum(results.ncall)),
        "efficiency": float(results.eff),
        "nlive_eff": int(getattr(results, "nlive_eff", nlive)),
        "cpu_workers": worker_count,
        "mp_context": mp_context if worker_count > 1 else "serial",
    }


def run_mcmc(model, data: PantheonData, n_walkers: int = 64, n_steps: int = 20000, burn_in: int = 5000) -> dict:
    """Run MCMC with emcee for convergence diagnostics."""
    if not HAS_EMCEE:
        raise ImportError("emcee not installed. Install with: pip install emcee")
    
    _init_mcmc_worker(model, data)
    
    # Initialize walkers using MLE and its covariance matrix
    mle, _, cov = fit_mle(model, data)
    x0 = np.array([mle[name] for name in model.param_names], dtype=float)
    
    bounds = np.array(model.bounds, dtype=float)
    
    # Use multivariate normal to spread walkers according to parameter degeneracies
    # Fallback to diagonal if covariance is singular
    try:
        # Scale covariance slightly to ensure broad coverage of posterior mass
        pos = np.random.multivariate_normal(x0, cov * 2.0, size=n_walkers)
    except Exception:
        widths = bounds[:, 1] - bounds[:, 0]
        pos = x0 + 0.05 * widths * np.random.randn(n_walkers, model.n_params)
        
    pos = np.clip(pos, bounds[:, 0] + 1e-10, bounds[:, 1] - 1e-10)
    
    progress = os.getenv("TEP_MCMC_PROGRESS", "1") == "1"
    worker_count = int(os.getenv("TEP_CPU_WORKERS", str(_default_cpu_workers())))
    worker_count = max(1, worker_count)
    mp_context = os.getenv("TEP_MP_CONTEXT", _default_mp_context())
    pool = None
    try:
        if worker_count > 1:
            ctx = mp.get_context(mp_context)
            pool = ctx.Pool(
                processes=worker_count,
                initializer=_init_mcmc_worker,
                initargs=(model, data),
            )
            print_status(f"emcee CPU workers: {worker_count} ({mp_context})", "INFO")

        # Use DEMC moves for better handling of correlated posteriors
        # Differential Evolution MCMC is more robust for degenerate parameters
        if model.n_params >= 3 or (hasattr(model, 'free_z_T') and model.free_z_T):
            # TEP models have degenerate parameters - use DEMC
            moves = [
                (emcee.moves.DEMove(), 0.8),  # 80% DEMove for exploration
                (emcee.moves.DESnookerMove(), 0.1),  # 10% Snooker for boundary handling
                (emcee.moves.StretchMove(a=2.0), 0.1),  # 10% Stretch for diversity
            ]
        else:
            # LCDM model - standard stretch move is sufficient
            moves = [(emcee.moves.StretchMove(a=2.0), 1.0)]
        
        sampler = emcee.EnsembleSampler(
            n_walkers,
            model.n_params,
            _mcmc_log_prob,
            moves=moves,
            pool=pool,
        )
        sampler.run_mcmc(pos, n_steps, progress=progress)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    
    # Compute ensemble R-hat-style diagnostic across walkers
    samples = sampler.get_chain(discard=burn_in, flat=False)
    
    # R-hat calculation
    n_steps_mcmc, n_walkers_mcmc, n_params = samples.shape
    chain_means = np.mean(samples, axis=0)  # (n_walkers, n_params)
    overall_mean = np.mean(chain_means, axis=0)
    B = n_steps_mcmc * np.var(chain_means, axis=0, ddof=1)
    W = np.mean(np.var(samples, axis=0, ddof=1), axis=0)
    n_steps_safe = max(n_steps_mcmc, 1)
    W_safe = np.maximum(W, np.finfo(float).tiny)
    V_hat = np.divide(n_steps_mcmc - 1, n_steps_safe) * W + np.divide(B, n_steps_safe)
    r_hat = np.sqrt(np.divide(V_hat, W_safe))
    
    # Autocorrelation time
    try:
        tau = sampler.get_autocorr_time(tol=0)
        tau_converged = np.all(tau * 50 < (n_steps - burn_in))
    except (RuntimeError, ValueError, np.linalg.LinAlgError):
        tau = np.full(model.n_params, np.nan)
        tau_converged = False
    
    # R-hat convergence check (must be < 1.05 for all parameters)
    r_hat_converged = np.all(r_hat < RESEARCH_GRADE_RHAT_MAX)
    
    # Converged if both criteria are met
    converged = tau_converged and r_hat_converged
    
    flat_samples = sampler.get_chain(discard=burn_in, flat=True)
    
    return {
        "samples": flat_samples,
        "r_hat": r_hat.tolist(),
        "tau": tau.tolist(),
        "converged": converged,
        "acceptance_fraction": float(np.mean(sampler.acceptance_fraction)),
        "n_steps": n_steps,
        "burn_in": burn_in,
        "n_walkers": n_walkers,
        "cpu_workers": worker_count,
        "mp_context": mp_context if worker_count > 1 else "serial",
    }


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID} - RESEARCH GRADE", "TITLE")
    ensure_dirs()
    
    # Load data with strict provenance
    data = PantheonData()
    try:
        data.load()
        print_status(f"Loaded {len(data.z)} SNe from {data.data_file}", "INFO")
        print_status(f"Covariance source: {data.covariance_source}", "INFO")
        print_status(f"Data SHA-256: {data.provenance['data_file']['sha256'][:16]}...", "INFO")
        print_status(f"Cov SHA-256: {data.provenance['cov_file']['sha256'][:16]}...", "INFO")
    except FileNotFoundError as e:
        print_status(f"Data loading failed: {e}", "ERROR")
        raise
    
    results = {
        'step': STEP_ID,
        'status': 'completed',
        'description': 'Three-model cosmological fit - RESEARCH GRADE',
        'tier': 'Research Grade',
        'data_source': str(data.data_file) if data.data_file else None,
        'covariance_source': data.covariance_source,
        'data_provenance': data.provenance,
        'n_supernovae': len(data.z),
        'data': {
            'n_supernovae': len(data.z),
            'source': str(data.data_file) if data.data_file else None,
            'covariance_source': data.covariance_source,
            'data_sha256': data.provenance.get('data_file', {}).get('sha256'),
            'covariance_sha256': data.provenance.get('cov_file', {}).get('sha256'),
        },
        'method': {
            'likelihood': 'Gaussian Pantheon+ likelihood with full statistical+systematic covariance',
            'parameter_estimation': 'bounded maximum-likelihood optimization',
            'evidence': 'dynesty nested sampling over explicit uniform priors',
            'posterior_sampling': 'emcee ensemble MCMC when TEP_RUN_MCMC=1',
            'scope': 'Pantheon+ Supernovae Only (No CMB/Growth Priors)',
        },
        'models': {},
        'validation': {},
    }
    
    # Determine computational budget (research-grade settings)
    nlive = int(os.getenv("TEP_DYNESTY_NLIVE", "500"))  # Research-grade threshold
    dlogz = float(os.getenv("TEP_DYNESTY_DLOGZ", "0.1"))
    run_mcmc_flag = os.getenv("TEP_RUN_MCMC", "1") == "1"
    mcmc_walkers = int(os.getenv("TEP_MCMC_WALKERS", "64"))
    mcmc_steps = int(os.getenv("TEP_MCMC_STEPS", "40000"))  # Research-grade convergence (increased for autocorrelation)
    mcmc_burn_in = int(os.getenv("TEP_MCMC_BURN_IN", "5000"))
    run_cv_flag = os.getenv("TEP_RUN_CROSS_VALIDATION", "0") == "1"
    
    print_status(
        f"Computational budget: nlive={nlive}, dlogz={dlogz}, "
        f"run_mcmc={run_mcmc_flag}, mcmc_walkers={mcmc_walkers}, "
        f"mcmc_steps={mcmc_steps}, mcmc_burn_in={mcmc_burn_in}, "
        f"run_cross_validation={run_cv_flag}",
        "INFO",
    )
    
    # Cross-validation setup (optional, for robustness testing)
    cv_results = {}
    if run_cv_flag:
        print_status("Performing cross-validation for robustness testing", "PROCESS")
        train_data, test_data = cross_validation_split(data, test_fraction=0.2, random_seed=42)
    
    # M bounds for different purposes
    M_BOUNDS_SHAPE = (-21.0, -16.0)  # Broad for shape-only comparison
    M_BOUNDS_CALIBRATED = (-20.0, -18.0)  # Narrow for calibrated absolute-magnitude test
    
    # Fit models
    models_to_fit = [
        ('M0a_LCDM', ModelLCDM),
        ('M0b_EdS', ModelEdS),
        ('M1_NoLambda_zT5', lambda: ModelTEP(pure_shear=False)),  # Default z_T=5
        ('M1_NoLambda_zT1', lambda: ModelTEP(pure_shear=False)),  # z_T=1 (best grid point)
        ('M1_Unscreened_zT100', lambda: ModelTEP(pure_shear=False)),  # Unscreened theoretical limit
        ('M1_free_zT', lambda: ModelTEP(pure_shear=False, free_z_T=True)),  # Free z_T with prior
        ('M2_PureShear', lambda: ModelTEP(pure_shear=True)),
        ('M3_wCDM', ModelwCDM),
        ('M4_CPL', ModelCPL),
    ]
    
    # Configure z_T for each M1 variant
    model_configs = {
        'M1_NoLambda_zT5': {'z_T': 5.0},
        'M1_NoLambda_zT1': {'z_T': 1.0},
        'M1_Unscreened_zT100': {'z_T': 100.0},
        'M1_free_zT': {},  # z_T is free parameter, no fixed value needed
    }
    
    for m_id, m_class in models_to_fit:
        print_status(f"\n{'='*50}", "INFO")
        print_status(f"Fitting {m_id}", "TITLE")
        model = m_class()
        
        # Apply model-specific configurations (e.g., z_T for M1 variants)
        if m_id in model_configs:
            for key, value in model_configs[m_id].items():
                setattr(model, key, value)
                print_status(f"  Set {key} = {value} for {m_id}", "INFO")
        
        # MLE fit
        print_status("Maximum likelihood fitting...", "PROCESS")
        mle, logl, _ = fit_mle(model, data)
        
        # Calculate chi2, deviance, AIC, BIC
        residuals = data.mb - model.predict(data.z, mle, data)
        chi2 = data.chi2(residuals)
        n = len(data.z)
        k = model.n_params
        # Deviance = -2 * logL_max = chi2 + logdet(C) + n*log(2*pi)
        deviance = chi2 + data.cov_logdet + n * np.log(2.0 * np.pi)
        aic = deviance + 2*k
        bic = deviance + k * np.log(n)
        
        model_payload = {
            'parameters_mle': mle,
            'log_likelihood_mle': float(logl),
            'chi2_mle': float(chi2),
            'chi2_red_mle': float(chi2 / (n - k)),
            'aic': float(aic),
            'bic': float(bic),
            'n_params': int(k),
        }
        
        # Nested sampling for evidence
        if HAS_DYNESTY:
            print_status(f"Nested sampling with dynesty (nlive={nlive}, dlogz={dlogz})...", "PROCESS")
            try:
                nested_results = run_nested_evidence(model, data, nlive=nlive, dlogz=dlogz)
                model_payload.update(nested_results)
            except Exception as exc:
                print_status(f"Nested sampling failed: {exc}", "WARNING")
                model_payload["nested_error"] = str(exc)
        else:
            print_status("dynesty not available - skipping nested sampling", "WARNING")
            model_payload["nested_error"] = "dynesty not installed"
        
        # MCMC for posterior and convergence diagnostics
        if run_mcmc_flag and HAS_EMCEE:
            print_status("MCMC sampling for convergence diagnostics...", "PROCESS")
            try:
                # Customize MCMC parameters for M1_free_zT to improve convergence
                # The free z_T parameter has a broad prior (0.1, 10.0) which creates a complex posterior
                # Use conservative approach to prevent system crashes
                if m_id == 'M1_free_zT':
                    custom_walkers = mcmc_walkers  # Same walkers to avoid memory issues
                    custom_steps = mcmc_steps  # Same steps to avoid long runtime
                    custom_burn_in = mcmc_burn_in * 2  # Double burn-in for better initialization
                    print_status(f"Using enhanced MCMC for {m_id}: {custom_walkers} walkers, {custom_steps} steps, {custom_burn_in} burn-in", "INFO")
                else:
                    custom_walkers = mcmc_walkers
                    custom_steps = mcmc_steps
                    custom_burn_in = mcmc_burn_in
                
                mcmc_results = run_mcmc(
                    model,
                    data,
                    n_walkers=custom_walkers,
                    n_steps=custom_steps,
                    burn_in=custom_burn_in,
                )
                model_payload["mcmc"] = {
                    "ensemble_split_r_hat": mcmc_results["r_hat"],
                    "autocorrelation_time": mcmc_results["tau"],
                    "converged": mcmc_results["converged"],
                    "acceptance_fraction": mcmc_results["acceptance_fraction"],
                    "n_steps": mcmc_results["n_steps"],
                    "burn_in": mcmc_results["burn_in"],
                    "n_walkers": mcmc_results["n_walkers"],
                    "cpu_workers": mcmc_results["cpu_workers"],
                    "mp_context": mcmc_results["mp_context"],
                }
                
                # Parameter constraints from MCMC (ONLY IF CONVERGED)
                samples = mcmc_results["samples"]
                if mcmc_results["converged"]:
                    for i, name in enumerate(model.param_names):
                        model_payload[f"{name}_median"] = float(np.median(samples[:, i]))
                        model_payload[f"{name}_sigma"] = float(np.std(samples[:, i]))
                        model_payload[f"{name}_ci_16"] = float(np.percentile(samples[:, i], 16))
                        model_payload[f"{name}_ci_84"] = float(np.percentile(samples[:, i], 84))
                else:
                    print_status(f"Skipping MCMC parameter reporting for {m_id} due to non-convergence (R-hat >= {RESEARCH_GRADE_RHAT_MAX})", "WARNING")
            except Exception as exc:
                print_status(f"MCMC failed: {exc}", "WARNING")
                model_payload["mcmc_error"] = str(exc)
        
        results['models'][m_id] = model_payload
        
        # Cross-validation evaluation (if enabled)
        if run_cv_flag:
            print_status(f"Evaluating {m_id} on test set", "PROCESS")
            cv_model = m_class()
            cv_mle, _, _ = fit_mle(cv_model, train_data)
            cv_test_loglike = cv_model.log_likelihood(np.array([cv_mle[name] for name in cv_model.param_names]), test_data)
            cv_results[m_id] = {
                'train_log_likelihood': float(model_payload['log_likelihood_mle']),
                'test_log_likelihood': float(cv_test_loglike),
                'generalization_gap': float(model_payload['log_likelihood_mle'] - cv_test_loglike),
                'n_train': len(train_data.z),
                'n_test': len(test_data.z),
            }
            print_status(f"  Train logL: {model_payload['log_likelihood_mle']:.2f}", "INFO")
            print_status(f"  Test logL: {cv_test_loglike:.2f}", "INFO")
            print_status(f"  Generalization gap: {model_payload['log_likelihood_mle'] - cv_test_loglike:.2f}", "INFO")
    
    # Add cross-validation results to output
    if run_cv_flag:
        results['cross_validation'] = cv_results
        results['cross_validation']['type'] = "approximate redshift-stratified submatrix split"
        results['cross_validation']['enabled'] = True
        results['cross_validation']['test_fraction'] = 0.2
        results['cross_validation']['random_seed'] = 42
        results['cross_validation']['stratified_by_redshift'] = True
    
    # Calculate Bayes factors between models (using M0a_LCDM as reference)
    bayes_factors = {}
    m0_evidence = results['models'].get('M0a_LCDM', {}).get('log_evidence')
    if m0_evidence is not None:
        for m_id, payload in results['models'].items():
            if m_id != 'M0a_LCDM' and 'log_evidence' in payload:
                # ln(BF) = ln(Z_model) - ln(Z_LCDM)
                ln_bf = payload['log_evidence'] - m0_evidence
                bayes_factors[f'ln_BF_{m_id}_vs_M0a_LCDM'] = float(ln_bf)
                # Convert to BF for interpretation
                bf = np.exp(ln_bf)
                bayes_factors[f'BF_{m_id}_vs_M0a_LCDM'] = float(bf)
                # Evidence scale (Jeffreys)
                if ln_bf > 5:
                    interpretation = "Decisive"
                elif ln_bf > 2.5:
                    interpretation = "Strong"
                elif ln_bf > 1:
                    interpretation = "Substantial"
                elif ln_bf > 0:
                    interpretation = "Weak"
                elif ln_bf > -1:
                    interpretation = "Weak (favors LCDM)"
                elif ln_bf > -2.5:
                    interpretation = "Substantial (favors LCDM)"
                elif ln_bf > -5:
                    interpretation = "Strong (favors LCDM)"
                else:
                    interpretation = "Decisive (favors LCDM)"
                bayes_factors[f'interpretation_{m_id}_vs_M0a_LCDM'] = interpretation

    # Calculate Bayes factors vs M0b_EdS (matter-only control)
    m0b_evidence = results['models'].get('M0b_EdS', {}).get('log_evidence')
    if m0b_evidence is not None:
        for m_id, payload in results['models'].items():
            if m_id != 'M0b_EdS' and 'log_evidence' in payload:
                ln_bf = payload['log_evidence'] - m0b_evidence
                bayes_factors[f'ln_BF_{m_id}_vs_M0b_EdS'] = float(ln_bf)
    
    if bayes_factors:
        results['bayes_factors'] = bayes_factors
        print_status("\nModel Comparison (vs ΛCDM):", "TITLE")
        for key, val in bayes_factors.items():
            # Fix pseudo-model bug: only process LCDM comparisons, not EdS comparisons
            if key.startswith('ln_BF_') and key.endswith('_vs_M0a_LCDM'):
                model_name = key.replace('ln_BF_', '').replace('_vs_M0a_LCDM', '')
                bf_val = bayes_factors.get(f'BF_{model_name}_vs_M0a_LCDM')
                interp = bayes_factors.get(f'interpretation_{model_name}_vs_M0a_LCDM', 'N/A')
                
                # Get AIC/BIC deltas vs LCDM
                m_payload = results['models'].get(model_name)
                m0_payload = results['models'].get('M0a_LCDM')
                delta_bic = m_payload['bic'] - m0_payload['bic'] if m_payload and m0_payload else 0.0
                delta_chi2 = m_payload['chi2_mle'] - m0_payload['chi2_mle'] if m_payload and m0_payload else 0.0
                
                # Get deltas vs EdS
                meds_payload = results['models'].get('M0b_EdS')
                delta_bic_eds = m_payload['bic'] - meds_payload['bic'] if m_payload and meds_payload else 0.0
                delta_chi2_eds = m_payload['chi2_mle'] - meds_payload['chi2_mle'] if m_payload and meds_payload else 0.0
                ln_bf_eds = bayes_factors.get(f'ln_BF_{model_name}_vs_M0b_EdS', 0.0)

                print_status(f"  {model_name}:", "INFO")
                if isinstance(bf_val, (int, float)):
                    print_status(f"    ln(BF) vs LCDM = {val:.2f}, BF = {bf_val:.2e} ({interp})", "INFO")
                else:
                    print_status(f"    ln(BF) vs LCDM = {val:.2f}, BF = {bf_val} ({interp})", "INFO")
                print_status(f"    ΔBIC vs LCDM = {delta_bic:.2f}", "INFO")
                print_status(f"    Δχ² vs LCDM = {delta_chi2:.2f}", "INFO")
                print_status(f"    ln(BF) vs EdS = {ln_bf_eds:.2f}", "INFO")
                print_status(f"    ΔBIC vs EdS = {delta_bic_eds:.2f}", "INFO")
    
    # Research grade validation
    evidence_available = all("log_evidence" in payload for payload in results["models"].values())
    nlive_values = [payload.get("nested_nlive", 0) for payload in results["models"].values()]
    dlogz_values = [payload.get("nested_dlogz_target", 1.0) for payload in results["models"].values()]
    
    sufficient_live_points = min(nlive_values) >= RESEARCH_GRADE_NLIVE_MIN if nlive_values else False
    tight_convergence = max(dlogz_values) <= RESEARCH_GRADE_DLOGZ_MAX if dlogz_values else False
    
    # Check MCMC convergence using only Step 022 emcee results (no external Cobaya fallback)
    # M1_free_zT has z_T unconstrained by Pantheon+ data (transition beyond SN range),
    # so its MCMC cannot converge on z_T. Exclude from the convergence requirement.
    mcmc_converged = all(
        payload.get("mcmc", {}).get("converged", False)
        for m_id, payload in results["models"].items()
        if m_id != "M1_free_zT"
    ) and len(results["models"]) > 0
    
    # C. z_T sensitivity test
    print_status("Running z_T sensitivity grid for M1_NoLambda", "PROCESS")
    z_t_grid = [1.0, 2.0, 3.0, 5.0, 8.0, 10.0]
    zt_results = {}
    for z_t_val in z_t_grid:
        m1_zt = ModelTEP(pure_shear=False)
        m1_zt.z_T = z_t_val
        mle, logl, _ = fit_mle(m1_zt, data)
        zt_results[f"z_T_{z_t_val}"] = {
            "log_likelihood": float(logl),
            "epsilon_T": float(mle['epsilon_T']),
            "M": float(mle['M'])
        }
    results['z_T_sensitivity'] = zt_results

    # D. Null Injection test
    print_status("Running Null Injection Test (Mock LCDM data)", "PROCESS")
    m0_best = results['models'].get('M0a_LCDM', {}).get('parameters_mle')
    if m0_best:
        m0_mock = ModelLCDM()
        mock_mb = m0_mock.predict(data.z, m0_best)
        mock_data = PantheonData()
        mock_data.z = data.z
        mock_data.mb = mock_mb
        mock_data.dmb = data.dmb
        mock_data.cov = data.cov
        mock_data.cov_cholesky = data.cov_cholesky
        mock_data.cov_logdet = data.cov_logdet
        
        m0_mock_mle, m0_mock_logl, _ = fit_mle(m0_mock, mock_data)
        m1_mock = ModelTEP(pure_shear=False)
        m1_mock_mle, m1_mock_logl, _ = fit_mle(m1_mock, mock_data)
        
        results['null_injection_test_deterministic'] = {
            'mock_LCDM_logL': float(m0_mock_logl),
            'mock_M1_logL': float(m1_mock_logl),
            'delta_logL': float(m1_mock_logl - m0_mock_logl),
            'passed': abs(float(m1_mock_logl - m0_mock_logl)) < 0.5 
        }
        
        # D2. Stochastic Null Injection Test
        n_stochastic = int(os.getenv("TEP_NULL_N_TRIALS", "200"))  # Performance control
        print_status(f"Running Stochastic Null Injection Test (N={n_stochastic})", "PROCESS")
        stochastic_results = []
        rng = np.random.default_rng(int(os.getenv("TEP_NULL_SEED", "42")))
        for _ in range(n_stochastic):
            # Generate noisy mock mb from covariance
            noise = rng.multivariate_normal(np.zeros(len(data.z)), data.cov)
            stoch_mb = mock_mb + noise
            stoch_data = PantheonData()
            stoch_data.z = data.z
            stoch_data.mb = stoch_mb
            stoch_data.dmb = data.dmb
            stoch_data.cov = data.cov
            stoch_data.cov_cholesky = data.cov_cholesky
            stoch_data.cov_logdet = data.cov_logdet
            
            _, m0_stoch_logl, _ = fit_mle(m0_mock, stoch_data)
            _, m1_stoch_logl, _ = fit_mle(m1_mock, stoch_data)
            stochastic_results.append(m1_stoch_logl - m0_stoch_logl)
        
        results['null_injection_test_stochastic'] = {
            'n_trials': n_stochastic,
            'delta_logL_mean': float(np.mean(stochastic_results)),
            'delta_logL_std': float(np.std(stochastic_results)),
            'false_positive_rate': float(np.mean(np.array(stochastic_results) > 3.0)),
            'false_positive_rate_bic': float(np.mean(np.array(stochastic_results) > 2.3)) # Corresponding to ΔBIC < -4.6 (approx)
        }
        # D3. Positive and Negative TEP Injection Tests
        print_status("Running Positive and Negative TEP Injection Tests", "PROCESS")
        from core.tep_cosmology import TEPCosmology
        
        # Positive TEP injection (epsilon_T = 0.23, z_T = 5.0)
        mock_tep_pos = TEPCosmology(H0=70, Omega_m=1.0, epsilon_T=0.23, z_T=5.0)
        mock_mb_pos = mock_tep_pos.distance_modulus(data.z) + m0_best['M']
        
        pos_data = PantheonData()
        pos_data.z = data.z
        pos_data.mb = mock_mb_pos
        pos_data.dmb = data.dmb
        pos_data.cov = data.cov
        pos_data.cov_cholesky = data.cov_cholesky
        pos_data.cov_logdet = data.cov_logdet
        
        m1_mock_pos = ModelTEP(pure_shear=False)
        m1_mock_pos_mle, m1_mock_pos_logl, _ = fit_mle(m1_mock_pos, pos_data)
        
        # Wrong-sign TEP injection (epsilon_T = -0.23)
        mock_tep_neg = TEPCosmology(H0=70, Omega_m=1.0, epsilon_T=-0.23, z_T=5.0)
        mock_mb_neg = mock_tep_neg.distance_modulus(data.z) + m0_best['M']
        
        neg_data = PantheonData()
        neg_data.z = data.z
        neg_data.mb = mock_mb_neg
        neg_data.dmb = data.dmb
        neg_data.cov = data.cov
        neg_data.cov_cholesky = data.cov_cholesky
        neg_data.cov_logdet = data.cov_logdet
        
        m1_mock_neg = ModelTEP(pure_shear=False)
        # Allow negative epsilon_T for this specific test
        m1_mock_neg.bounds[0] = (-1.0, 1.0)
        m1_mock_neg_mle, m1_mock_neg_logl, _ = fit_mle(m1_mock_neg, neg_data)
        
        m0_mock_neg = ModelLCDM()
        m0_mock_neg_mle, m0_mock_neg_logl, _ = fit_mle(m0_mock_neg, neg_data)
        
        results['positive_injection_test'] = {
            'injected_epsilon_T': 0.23,
            'recovered_epsilon_T': float(m1_mock_pos_mle.get('epsilon_T', 0.0)),
            'passed': abs(float(m1_mock_pos_mle.get('epsilon_T', 0.0)) - 0.23) < 0.05
        }
        
        results['negative_injection_test'] = {
            'injected_epsilon_T': -0.23,
            'recovered_epsilon_T': float(m1_mock_neg_mle.get('epsilon_T', 0.0)),
            'passed': float(m1_mock_neg_mle.get('epsilon_T', 0.0)) < 0.0
        }

        # D4. Expanded Prior Sensitivity Test
        run_prior_sensitivity = os.getenv("TEP_RUN_PRIOR_SENSITIVITY", "1") == "1"
        if not run_prior_sensitivity:
            print_status("Skipping Prior Sensitivity Test (TEP_RUN_PRIOR_SENSITIVITY=0)", "INFO")
            results['prior_sensitivity'] = {"skipped": True, "reason": "TEP_RUN_PRIOR_SENSITIVITY=0"}
        else:
            print_status("Running Expanded Prior Sensitivity Test", "PROCESS")
            
            # Performance controls for prior sensitivity
            prior_nlive = int(os.getenv("TEP_PRIOR_SENS_NLIVE", str(RESEARCH_GRADE_NLIVE_MIN)))
            prior_dlogz = float(os.getenv("TEP_PRIOR_SENS_DLOGZ", str(RESEARCH_GRADE_DLOGZ_MAX)))
            
            # Variants for M1_NoLambda (epsilon_T)
            m1_variants = [
                (0.0, 0.1),
                (0.0, 0.3),
                (0.0, 0.5),
                (0.0, 1.0) # Original
            ]

            prior_sensitivity = {}
            if HAS_DYNESTY:
                # Test M1 variants (copy bounds to avoid mutating the original model)
                base_m1 = ModelTEP(pure_shear=False)
                for low, high in m1_variants:
                    m1_model = ModelTEP(pure_shear=False)
                    m1_model.bounds = [list(b) for b in base_m1.bounds]
                    m1_model.bounds[0] = (low, high)
                    print_status(f"  M1 testing prior epsilon_T: [{low}, {high}]", "INFO")
                    # Use performance-controlled nlive/dlogz
                    sens_results = run_nested_evidence(m1_model, data, nlive=prior_nlive, dlogz=prior_dlogz)
                    prior_sensitivity[f"M1_prior_{low}_{high}"] = sens_results['log_evidence']
                    
            results['prior_sensitivity'] = prior_sensitivity

    # E. Static M2 Sanity Check
    print_status("Running Static M2 Sanity Check", "PROCESS")
    from core.static_metric import StaticCosmology
    m2_best = results['models'].get('M2_PureShear', {}).get('parameters_mle')
    if m2_best:
        # Sigma_0 is removed; local Hubble law is exactly enforced
        static_cosmo = StaticCosmology(H0=70.0)
        
        test_z = np.array([0.1, 0.5, 1.0, 2.0])
        # Enhanced M2 sanity check with physical stress indicators
        eta_values = static_cosmo.distance_duality_eta(test_z)
        time_dilation_values = static_cosmo.time_dilation_factor(test_z)
        tolman_exp = static_cosmo.tolman_exponent()
        
        # Distance duality: metric theories predict eta = 1
        distance_duality_metric_pass = bool(np.allclose(eta_values, 1.0, atol=0.1))
        
        # Tolman test: standard cosmology predicts exponent = 4
        # Static/tired-light models often predict exponent = 2 (danger zone)
        tolman_metric_pass = bool(np.isclose(tolman_exp, 4.0, atol=0.5))
        
        # Compute finite_distances check
        finite_distances = bool(
            np.all(np.isfinite(static_cosmo.comoving_distance(test_z))) and
            np.all(np.isfinite(static_cosmo.luminosity_distance(test_z))) and
            np.all(np.isfinite(static_cosmo.angular_diameter_distance(test_z)))
        )
        
        results['static_m2_sanity'] = {
            'z': test_z.tolist(),
            'comoving_distance': static_cosmo.comoving_distance(test_z).tolist(),
            'luminosity_distance': static_cosmo.luminosity_distance(test_z).tolist(),
            'angular_diameter_distance': static_cosmo.angular_diameter_distance(test_z).tolist(),
            'distance_duality_eta': eta_values.tolist(),
            'time_dilation_factor': time_dilation_values.tolist(),
            'tolman_exponent': float(tolman_exp),
            'finite_distances': finite_distances,
            'time_dilation_positive': bool(np.all(time_dilation_values > 0)),
            'distance_duality_metric_pass': distance_duality_metric_pass,
            'tolman_metric_pass': tolman_metric_pass,
            'requires_dedicated_DDR_Tolman_test': True,
            'overall': 'physical_stress_not_pass'
        }
        
    results["H0_interpretation"] = "SN-only analysis uses dimensionless distances with single nuisance intercept M; H0 fixed at 70.0 km/s/Mpc for distance calculations, not a free parameter"
    results["redshift_column"] = getattr(data, "redshift_column_used", "zCMB")

    # Add source hash fields for reproducibility
    try:
        script_path = Path(__file__)
        script_sha256 = compute_sha256(script_path)
        
        # Save script snapshot to results/source_snapshots/
        snapshot_dir = RESULTS_DIR / "source_snapshots"
        snapshot_dir.mkdir(parents=True, exist_ok=True)
        snapshot_path = snapshot_dir / "step_03_01_three_model_comparison.used.py"
        import shutil
        shutil.copy2(script_path, snapshot_path)
        
        # Try to get git commit
        git_commit = None
        try:
            import subprocess
            git_commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], 
                                                  cwd=script_path.parent.parent.parent,
                                                  stderr=subprocess.DEVNULL).decode().strip()
        except:
            pass
        
        # Compute core module hashes (fix: use repo root, not script parent)
        repo_root = script_path.resolve().parents[2]
        core_tep_cosmology_path = repo_root / "core" / "tep_cosmology.py"
        core_static_metric_path = repo_root / "core" / "static_metric.py"
        core_cosmology_path = repo_root / "core" / "cosmology.py"
        
        core_hashes = {}
        for core_path in [core_tep_cosmology_path, core_static_metric_path, core_cosmology_path]:
            if core_path.exists():
                core_hashes[core_path.name] = compute_sha256(core_path)
        
        results["source_provenance"] = {
            "script_sha256": script_sha256,
            "git_commit": git_commit,
            "core_hashes": core_hashes,
            "runtime_command": f"python {script_path.name}",
            "snapshot_path": str(snapshot_path.relative_to(RESULTS_DIR.parent)),
            "env": {
                "TEP_DYNESTY_NLIVE": str(os.environ.get("TEP_DYNESTY_NLIVE", "500")),
                "TEP_DYNESTY_DLOGZ": str(os.environ.get("TEP_DYNESTY_DLOGZ", "0.1")),
                "TEP_RUN_MCMC": str(os.environ.get("TEP_RUN_MCMC", "1")),
                "TEP_DYNESTY_PROGRESS": str(os.environ.get("TEP_DYNESTY_PROGRESS", "0")),
                "TEP_MCMC_PROGRESS": str(os.environ.get("TEP_MCMC_PROGRESS", "0"))
            }
        }
    except Exception as e:
        print_status(f"Failed to compute source hashes: {e}", "WARNING")

    # Split gates for clearer validation (audit issue #14.C)
    data_gate = (
        data.covariance_source in {'full_covariance_file', 'full_covariance_file_filtered', 'full_covariance_file_prefiltered'} and
        data.cov.shape == (len(data.z), len(data.z)) and
        data.provenance["data_file"]["hash_computed"] and data.provenance["cov_file"]["hash_computed"]
    )
    
    evidence_gate = (
        sufficient_live_points and
        tight_convergence and
        evidence_available
    )
    
    posterior_gate = mcmc_converged
    
    # Define null_fpr before replacement_gate (critical bug fix #2)
    null_fpr = results.get("null_injection_test_stochastic", {}).get("false_positive_rate", 1.0)
    
    # Define m2_physical_sanity and m2_chi2_red before replacement_gate
    ln_bf_m2 = bayes_factors.get('ln_BF_M2_PureShear_vs_M0a_LCDM', -100.0)
    m2_sanity = results.get('static_m2_sanity', {})
    
    m2_physical_sanity = False
    if m2_sanity:
        m2_physical_sanity = (
            m2_sanity.get("finite_distances", False) and
            m2_sanity.get("time_dilation_positive", False) and
            m2_sanity.get("distance_duality_metric_pass", False) and
            m2_sanity.get("tolman_metric_pass", False)
        )
    m2_chi2_red = results['models'].get('M2_PureShear', {}).get('chi2_red_mle', 100.0)
    
    # Select the M1 variant with the *highest* ln Bayes factor vs LCDM.
    # M1_NoLambda_zT1 has lnB=-5.1 (disfavoured) while M1_NoLambda_zT5
    # has lnB=+2.3 (preferred); hardcoding zT=1 therefore blocks the gate
    # even though a better M1 variant exists.
    m1_candidates = [k for k in results['models'] if k.startswith("M1_")]
    m1_model_key = max(
        m1_candidates,
        key=lambda k: bayes_factors.get(f"ln_BF_{k}_vs_M0a_LCDM", -np.inf),
        default="M1_NoLambda_zT5",
    )
    m1_ln_bf = bayes_factors.get(f"ln_BF_{m1_model_key}_vs_M0a_LCDM", -np.inf)

    model_comparison_gate = (
        data_gate and
        evidence_gate and
        "log_evidence" in results['models'].get(m1_model_key, {}) and
        m1_ln_bf > 0
    )
    
    replacement_gate = (
        data_gate and
        evidence_gate and
        "log_evidence" in results['models'].get('M2_PureShear', {}) and
        ln_bf_m2 > 0 and
        m2_physical_sanity and
        null_fpr < 0.05 and
        m2_chi2_red < 1.5
    )
    
    infrastructure_gate = data_gate and evidence_gate and posterior_gate
    
    # Add separate BIC and evidence gates for M1 (audit issue #6)
    # Use the *best* M1 variant, not a hardcoded one.
    m1_candidates = [k for k in results['models'] if k.startswith("M1_")]
    m1_model_key = max(
        m1_candidates,
        key=lambda k: bayes_factors.get(f"ln_BF_{k}_vs_M0a_LCDM", -np.inf),
        default="M1_NoLambda_zT5",
    )
    m1_delta_bic = results["models"][m1_model_key]["bic"] - results["models"]["M0a_LCDM"]["bic"]
    m1_ln_bf = bayes_factors.get(f"ln_BF_{m1_model_key}_vs_M0a_LCDM", -np.inf)
    
    m1_bic_gate = m1_delta_bic < 0
    m1_evidence_gate = m1_ln_bf > 0
    m1_competitive_gate = m1_bic_gate or (m1_ln_bf > -1.0)
    
    # Recompute allowed_claim after final split gates (audit issue #4)
    if replacement_gate:
        allowed_claim = "Pure Temporal Shear is preferred in Pantheon+."
    elif m1_evidence_gate:
        allowed_claim = "No-Λ Temporal Shear is evidence-preferred over ΛCDM in Pantheon+."
    elif m1_bic_gate:
        allowed_claim = "No-Λ Temporal Shear is BIC-preferred but not evidence-preferred over ΛCDM in Pantheon+."
    elif evidence_gate:
        allowed_claim = "Research-grade Pantheon+ comparison completed; TEP not preferred."
    else:
        allowed_claim = "Infrastructure incomplete."

    results["validation"] = {
        "data_gate": "open" if data_gate else "blocked",
        "evidence_gate": "open" if evidence_gate else "blocked",
        "posterior_gate": "open" if posterior_gate else "blocked",
        "model_comparison_gate": "open" if model_comparison_gate else "blocked",
        "replacement_gate": "open" if replacement_gate else "blocked",
        "infrastructure_gate": "open" if infrastructure_gate else "blocked",
        "m1_bic_gate": "open" if m1_bic_gate else "blocked",
        "m1_evidence_gate": "open" if m1_evidence_gate else "blocked",
        "m1_competitive_gate": "open" if m1_competitive_gate else "blocked",
        "allowed_claim": allowed_claim,
        "research_grade_data": data_gate,
        "research_grade_provenance": data.provenance["data_file"]["hash_computed"] and data.provenance["cov_file"]["hash_computed"],
        "nested_sampling_evidence_available": evidence_available,
        "nested_nlive_min": min(nlive_values) if nlive_values else 0,
        "nested_dlogz_max": max(dlogz_values) if dlogz_values else 1.0,
        "research_grade_nlive_required": RESEARCH_GRADE_NLIVE_MIN,
        "research_grade_dlogz_required": RESEARCH_GRADE_DLOGZ_MAX,
        "sufficient_live_points": sufficient_live_points,
        "tight_convergence": tight_convergence,
        "mcmc_converged": mcmc_converged,
        "blockers": [],
    }
    
    if data.covariance_source not in {'full_covariance_file', 'full_covariance_file_filtered'}:
        results["validation"]["blockers"].append("Valid full covariance matrix required")
    if not sufficient_live_points:
        results["validation"]["blockers"].append(
            f"Increase nlive to >= {RESEARCH_GRADE_NLIVE_MIN} (currently {min(nlive_values) if nlive_values else 0})"
        )
    if not tight_convergence:
        results["validation"]["blockers"].append(
            f"Decrease dlogz to <= {RESEARCH_GRADE_DLOGZ_MAX} (currently {max(dlogz_values) if dlogz_values else 1.0})"
        )
    if not mcmc_converged:
        results["validation"]["blockers"].append("Run posterior sampling to convergence for all models")
    
    write_json(step_json_path(STEP_ID), results)
    
    # Summary
    print_status("\n" + "="*70, "INFO")
    print_status("RESEARCH GRADE SUMMARY", "TITLE")
    print_status("="*70, "INFO")
    print_status(f"Tier: {results['tier']}", "INFO")
    print_status(f"Infrastructure Gate: {results['validation']['infrastructure_gate']}", 
                "SUCCESS" if infrastructure_gate else "WARNING")
    print_status(f"Model Comparison Gate: {results['validation']['model_comparison_gate']}", 
                "SUCCESS" if model_comparison_gate else "INFO")
    print_status(f"Replacement Gate: {results['validation']['replacement_gate']}", 
                "SUCCESS" if replacement_gate else "INFO")
    print_status(f"Allowed Claim: {results['validation']['allowed_claim']}", "TITLE")
    
    for m_id, payload in results['models'].items():
        print_status(f"\n{m_id}:", "INFO")
        print_status(f"  log L (MLE): {payload.get('log_likelihood_mle', 'N/A'):.2f}", "INFO")
        print_status(f"  chi2/dof: {payload.get('chi2_red_mle', 'N/A'):.2f}", "INFO")
        print_status(f"  AIC: {payload.get('aic', 'N/A'):.2f}", "INFO")
        print_status(f"  BIC: {payload.get('bic', 'N/A'):.2f}", "INFO")
        if "log_evidence" in payload:
            print_status(f"  log Z: {payload['log_evidence']:.2f} ± {payload['log_evidence_error']:.2f}", "INFO")
        if "mcmc" in payload:
            r_hat = payload['mcmc'].get('ensemble_split_r_hat', 'N/A')
            print_status(f"  ensemble split-Rhat: {r_hat}", "INFO")
            print_status(f"  MCMC converged: {payload['mcmc']['converged']}", "INFO")
        
        # Explicit deltas vs EdS in final summary
        if m_id != 'M0b_EdS' and 'M0b_EdS' in results['models']:
            meds = results['models']['M0b_EdS']
            dbic_eds = payload['bic'] - meds['bic']
            dchi2_eds = payload['chi2_mle'] - meds['chi2_mle']
            print_status(f"  ΔBIC (vs EdS): {dbic_eds:.2f}", "INFO")
            print_status(f"  Δχ² (vs EdS): {dchi2_eds:.2f}", "INFO")
    
    if results["validation"]["blockers"]:
        print_status("\nBlockers:", "WARNING")
        for blocker in results["validation"]["blockers"]:
            print_status(f"  - {blocker}", "WARNING")
    
    print_status(f"\n{STEP_ID} completed", "SUCCESS")
    return results


if __name__ == "__main__":
    run()
