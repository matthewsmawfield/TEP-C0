# TEP-C0 Logging System

## Overview

The TEP-C0 pipeline provides two levels of logging:

1. **Standard Logging** (`logger.py`) - Production logging with color-coded console output
2. **Verbose Logging** (`verbose_logger.py`) - Highly detailed logging for debugging and analysis

## Standard Logging (`logger.py`)

### Features
- Color-coded console output (green=success, red=error, blue=process, etc.)
- File logging without color codes
- Custom log levels: PROCESS (25), SUCCESS (26), TEST (27)
- TEPFormatter with ANSI color support
- TEPFileFormatter for clean file output

### Usage
```python
from scripts.utils.logger import TEPLogger, print_status

logger = TEPLogger("my_step", log_file_path=Path("logs/my_step.log"))
print_status("Starting operation", "PROCESS")
print_status("Operation complete", "SUCCESS")
```

### Log Levels
| Level | Color | Purpose |
|-------|-------|---------|
| INFO | White | General information |
| PROCESS | Blue | Ongoing operations |
| SUCCESS | Green | Completed operations |
| WARNING | Yellow | Non-fatal issues |
| ERROR | Red | Errors |
| TEST | Magenta | Test results |
| DEBUG | Gray | Debug information |
| TITLE | Cyan | Section titles |

## Verbose Logging (`verbose_logger.py`)

### Features
- **Operation tracking**: Every function call logged with timing
- **Parameter logging**: All parameters tracked with section hierarchy
- **Data summaries**: Automatic summary statistics for arrays
- **Performance metrics**: Execution time in milliseconds
- **Environment capture**: Python version, platform, paths
- **Structured JSON output**: Machine-readable execution record
- **Error tracebacks**: Full stack traces on exceptions

### Output Files
| File | Purpose |
|------|---------|
| `logs/verbose/{step_id}_verbose.log` | Detailed human-readable log |
| `results/outputs/verbose/{step_id}_verbose.json` | Structured machine-readable record |

### Usage

#### Option 1: Context Manager
```python
from scripts.utils.verbose_logger import log_step_verbose

with log_step_verbose("step_033", "Cobaya TEP Inference") as logger:
    logger.log_checkpoint("Loading data")
    data = logger.log_operation("Load Pantheon+", load_pantheon_data)
    logger.log_param("n_sne", len(data.z))
    logger.log_data_summary("redshifts", data.z)
```

#### Option 2: Direct Usage
```python
from scripts.utils.verbose_logger import VerboseStepLogger

logger = VerboseStepLogger("step_033", verbose=True)
logger.start_step("Cobaya TEP Inference")

try:
    data = logger.log_operation("Load data", load_func)
    logger.log_param("n_sne", 1701)
    logger.log_metric("load_time", 45.2, "ms")
    logger.end_step("success")
except Exception as e:
    logger.end_step("error", str(e))
```

### JSON Output Structure
```json
{
  "step_id": "step_033",
  "start_time": "2026-05-02T12:00:00",
  "end_time": "2026-05-02T12:05:32",
  "duration_ms": 32000,
  "status": "success",
  "operations": [
    {
      "name": "Load Pantheon+",
      "duration_ms": 150,
      "status": "success",
      "params": {
        "file": "data/raw/Pantheon+SH0ES.dat"
      }
    }
  ],
  "parameters": {
    "n_sne": 1701,
    "redshift_range": [0.0012, 2.26]
  },
  "environment": {
    "python_version": "3.13.0",
    "platform": "darwin",
    "cpu_count": "8"
  }
}
```

## Integration with Pipeline Steps

### Example: Enhanced Step with Verbose Logging
```python
def run() -> dict:
    if HAS_VERBOSE_LOGGER:
        with log_step_verbose("step_XXX") as vlogger:
            return _run_with_logger(vlogger)
    else:
        return _run_basic()

def _run_with_logger(vlogger):
    vlogger.start_step("Step Description")
    
    # Check dependencies
    result = vlogger.log_operation("Check CLASS", check_tep_class)
    vlogger.log_param("class.available", result[0])
    
    # Load data
    data = vlogger.log_operation("Load data", load_pantheon_data)
    vlogger.log_data_summary("z", data.z)
    
    # Run MCMC
    vlogger.log_checkpoint("Starting MCMC")
    mcmc = vlogger.log_operation("Run MCMC", run_mcmc, data)
    vlogger.log_metric("mcmc_time", 120000, "ms")
    
    vlogger.end_step("success")
    return {"status": "success"}
```

## Environment Variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `TEP_VERBOSE` | `1` | Enable verbose logging globally |
| `TEP_COBAYA_SAMPLES` | `5000` | Max samples for Cobaya |
| `TEP_COBAYA_TIMEOUT` | `300` | Timeout in seconds |

## Migration Guide

### From Basic Print
```python
# Before
print(f"Loaded {n} SNe")

# After (verbose)
logger.log_param("n_sne", n)
logger.info(f"Loaded {n} SNe")
```

### From Standard Logger
```python
# Before
print_status("Loading data", "PROCESS")
data = load_data()
print_status(f"Loaded {len(data)} items", "SUCCESS")

# After (verbose)
vlogger.log_checkpoint("Loading data")
data = vlogger.log_operation("Load data", load_data)
vlogger.log_data_summary("data", data)
```

## File Locations

```
logs/
├── verbose/                        # Verbose logs
│   ├── step_033_cobaya_tep_inference_verbose.log
│   └── step_033_cobaya_tep_inference_verbose.json
└── step_XXX.log                   # Standard logs

results/outputs/
├── verbose/                        # Verbose JSON outputs
│   └── step_033_cobaya_tep_inference_verbose.json
└── step_XXX.json                  # Standard JSON outputs
```

## Best Practices

1. **Always use verbose logging in development** for detailed debugging
2. **Use standard logging in production** for cleaner output
3. **Wrap operations in `log_operation()`** for automatic timing
4. **Log parameters with `log_param()`** for configuration tracking
5. **Use `log_checkpoint()`** for major milestones
6. **Save JSON outputs** for programmatic analysis
7. **Set `verbose=False`** for standard runs to reduce output

## Testing

```python
# Test verbose logger
python -c "from scripts.utils.verbose_logger import VerboseStepLogger; print('OK')"

# Test standard logger
python -c "from scripts.utils.logger import TEPLogger; print('OK')"

# Run verbose step
python scripts/steps/step_033_cobaya_tep_inference_verbose.py
```

## Date

2026-05-02 - Verbose logging system implemented
