#!/usr/bin/env python3
"""
TEP-C0 Highly Verbose Logging System
=====================================

Provides comprehensive logging for every step with:
- Timestamped detailed operation logs
- Parameter tracking at every stage
- Performance metrics (timing, memory)
- Debug information for every function call
- Automatic output file generation

Usage:
    from verbose_logger import VerboseStepLogger
    
    logger = VerboseStepLogger("step_033", verbose=True)
    logger.start_step()
    logger.log_operation("Loading data", lambda: load_data())
    logger.log_param("n_samples", 1000)
    logger.end_step()
"""

from __future__ import annotations

import json
import logging
import sys
import time
import traceback
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union
import os

# Import base logger
sys.path.insert(0, str(Path(__file__).parent))
from logger import TEPLogger, TEPFormatter, TEPFileFormatter, print_status

PACKAGE_ROOT = Path(__file__).resolve().parents[2]


@dataclass
class OperationRecord:
    """Record of a single operation within a step."""
    name: str
    start_time: float
    end_time: Optional[float] = None
    status: str = "running"
    error: Optional[str] = None
    params: Dict[str, Any] = field(default_factory=dict)
    result_summary: Optional[str] = None
    
    @property
    def duration_ms(self) -> float:
        if self.end_time is None:
            return 0.0
        return (self.end_time - self.start_time) * 1000


@dataclass
class StepRecord:
    """Complete record of a step execution."""
    step_id: str
    start_time: float
    end_time: Optional[float] = None
    status: str = "running"
    operations: List[OperationRecord] = field(default_factory=list)
    parameters: Dict[str, Any] = field(default_factory=dict)
    environment: Dict[str, str] = field(default_factory=dict)
    error: Optional[str] = None
    
    @property
    def duration_ms(self) -> float:
        if self.end_time is None:
            return (time.time() - self.start_time) * 1000
        return (self.end_time - self.start_time) * 1000
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "step_id": self.step_id,
            "start_time": datetime.fromtimestamp(self.start_time).isoformat(),
            "end_time": datetime.fromtimestamp(self.end_time).isoformat() if self.end_time else None,
            "duration_ms": self.duration_ms,
            "status": self.status,
            "operations": [
                {
                    "name": op.name,
                    "duration_ms": op.duration_ms,
                    "status": op.status,
                    "error": op.error,
                    "params": op.params,
                    "result_summary": op.result_summary,
                }
                for op in self.operations
            ],
            "parameters": self.parameters,
            "environment": self.environment,
            "error": self.error,
        }


class VerboseStepLogger:
    """Highly verbose logger for pipeline steps.
    
    Provides detailed logging for every operation, automatic timing,
    parameter tracking, and comprehensive output generation.
    
    Attributes:
        step_id: Unique identifier for this step
        verbose: Whether to enable verbose logging
        log_file: Path to detailed log file
        json_output: Path to structured JSON output
    """
    
    def __init__(
        self,
        step_id: str,
        verbose: bool = True,
        log_dir: Optional[Path] = None,
        output_dir: Optional[Path] = None,
    ):
        self.step_id = step_id
        self.verbose = verbose
        self.record = StepRecord(
            step_id=step_id,
            start_time=time.time(),
            environment=self._collect_environment(),
        )
        
        # Setup paths
        self.log_dir = log_dir or (PACKAGE_ROOT / "logs" / "verbose")
        self.output_dir = output_dir or (PACKAGE_ROOT / "results" / "outputs" / "verbose")
        
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.log_file = self.log_dir / f"{step_id}_verbose.log"
        self.json_output = self.output_dir / f"{step_id}_verbose.json"
        
        # Setup logging
        self._setup_logging()
        
        # Current operation tracking
        self._current_operation: Optional[OperationRecord] = None
    
    def _collect_environment(self) -> Dict[str, str]:
        """Collect relevant environment information."""
        return {
            "timestamp": datetime.now().isoformat(),
            "python_version": sys.version,
            "platform": sys.platform,
            "cwd": str(Path.cwd()),
            "tep_class_path": os.environ.get("TEP_CLASS_PYTHONPATH", "not_set"),
            "cpu_count": str(os.cpu_count() or "unknown"),
        }
    
    def _setup_logging(self):
        """Setup verbose logging to file."""
        self.logger = logging.getLogger(f"verbose_{self.step_id}")
        self.logger.setLevel(logging.DEBUG if self.verbose else logging.INFO)
        self.logger.handlers.clear()
        
        # Console handler with TEP formatting
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG if self.verbose else logging.INFO)
        ch.setFormatter(TEPFormatter())
        self.logger.addHandler(ch)
        
        # File handler for detailed log
        fh = logging.FileHandler(self.log_file, mode='w', encoding='utf-8')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(TEPFileFormatter())
        self.logger.addHandler(fh)
    
    def start_step(self, description: Optional[str] = None):
        """Begin step logging."""
        desc = description or f"Step {self.step_id}"
        self.logger.info("")
        self.logger.info("=" * 80)
        self.logger.info(f"STARTING: {desc}")
        self.logger.info("=" * 80)
        self.logger.info(f"Log file: {self.log_file}")
        self.logger.info(f"JSON output: {self.json_output}")
        self.logger.info(f"Environment: {json.dumps(self.record.environment, indent=2)}")
        self.logger.info("")
        
        # Set global step logger
        try:
            from logger import set_step_logger
            base_logger = TEPLogger(self.step_id, log_file_path=self.log_file)
            set_step_logger(base_logger)
        except ImportError:
            pass
    
    def log_param(self, name: str, value: Any, section: Optional[str] = None):
        """Log a parameter value."""
        key = f"{section}.{name}" if section else name
        
        # Convert value to serializable format
        if isinstance(value, (int, float, str, bool)):
            serializable_value = value
        elif isinstance(value, (list, tuple)):
            serializable_value = list(value)[:10]  # Truncate long lists
            if len(value) > 10:
                serializable_value.append(f"... ({len(value)} total)")
        elif isinstance(value, dict):
            serializable_value = {k: str(v)[:100] for k, v in list(value.items())[:10]}
        else:
            serializable_value = str(value)[:100]
        
        self.record.parameters[key] = serializable_value
        
        if self.verbose:
            self.logger.debug(f"[PARAM] {key} = {serializable_value}")
    
    def log_operation(self, name: str, func: Callable, *args, **kwargs) -> Any:
        """Execute and log an operation with timing."""
        op = OperationRecord(name=name, start_time=time.time())
        self._current_operation = op
        self.record.operations.append(op)
        
        self.logger.info(f"[OPERATION] Starting: {name}")
        
        try:
            # Log input parameters
            if args:
                op.params["args_count"] = len(args)
                op.params["args_types"] = [type(a).__name__ for a in args]
            if kwargs:
                op.params["kwargs_keys"] = list(kwargs.keys())
            
            # Execute function
            result = func(*args, **kwargs)
            
            # Record success
            op.end_time = time.time()
            op.status = "success"
            op.result_summary = self._summarize_result(result)
            
            self.logger.success(f"[OPERATION] Completed: {name} ({op.duration_ms:.1f}ms)")
            
            return result
            
        except Exception as e:
            op.end_time = time.time()
            op.status = "error"
            op.error = str(e)
            op.params["traceback"] = traceback.format_exc()
            
            self.logger.error(f"[OPERATION] Failed: {name} - {e}")
            raise
        
        finally:
            self._current_operation = None
    
    def _summarize_result(self, result: Any) -> str:
        """Create a summary of a result."""
        if result is None:
            return "None"
        elif isinstance(result, (int, float, str, bool)):
            return str(result)
        elif isinstance(result, (list, tuple)):
            return f"{type(result).__name__} with {len(result)} items"
        elif isinstance(result, dict):
            return f"dict with {len(result)} keys"
        elif hasattr(result, '__class__'):
            return f"{result.__class__.__name__} instance"
        else:
            return str(type(result).__name__)
    
    def log_data_summary(self, name: str, data: Any):
        """Log summary statistics for data."""
        summary = {}
        
        if hasattr(data, '__len__'):
            summary['length'] = len(data)
        
        if hasattr(data, 'shape'):
            summary['shape'] = data.shape
        
        if hasattr(data, 'dtype'):
            summary['dtype'] = str(data.dtype)
        
        if hasattr(data, 'min') and hasattr(data, 'max'):
            try:
                summary['min'] = float(data.min())
                summary['max'] = float(data.max())
            except (TypeError, ValueError):
                pass
        
        if hasattr(data, 'mean'):
            try:
                summary['mean'] = float(data.mean())
            except (TypeError, ValueError):
                pass
        
        self.log_param(f"data_summary.{name}", summary)
        self.logger.debug(f"[DATA] {name}: {summary}")
    
    def log_checkpoint(self, message: str):
        """Log a checkpoint message."""
        self.logger.info(f"[CHECKPOINT] {message}")
    
    def log_metric(self, name: str, value: float, unit: str = ""):
        """Log a performance metric."""
        unit_str = f" {unit}" if unit else ""
        self.log_param(f"metric.{name}", f"{value}{unit_str}")
        self.logger.info(f"[METRIC] {name}: {value}{unit_str}")
    
    def end_step(self, status: str = "success", error: Optional[str] = None):
        """End step logging and save outputs."""
        self.record.end_time = time.time()
        self.record.status = status
        self.record.error = error
        
        # Log completion
        self.logger.info("")
        self.logger.info("=" * 80)
        self.logger.info(f"COMPLETED: {self.step_id}")
        self.logger.info(f"Status: {status.upper()}")
        self.logger.info(f"Duration: {self.record.duration_ms:.1f}ms")
        self.logger.info(f"Operations: {len(self.record.operations)}")
        if error:
            self.logger.error(f"Error: {error}")
        self.logger.info("=" * 80)
        
        # Save JSON output
        self._save_json_output()
        
        return self.record.to_dict()
    
    def _save_json_output(self):
        """Save verbose output to JSON."""
        try:
            output_data = self.record.to_dict()
            self.json_output.write_text(
                json.dumps(output_data, indent=2, default=str), 
                encoding='utf-8'
            )
            self.logger.info(f"[OUTPUT] Saved verbose JSON to {self.json_output}")
        except Exception as e:
            self.logger.error(f"[OUTPUT] Failed to save JSON: {e}")
    
    def info(self, message: str):
        """Log info message."""
        self.logger.info(message)
    
    def debug(self, message: str):
        """Log debug message (only if verbose)."""
        self.logger.debug(message)
    
    def warning(self, message: str):
        """Log warning message."""
        self.logger.warning(message)
    
    def error(self, message: str):
        """Log error message."""
        self.logger.error(message)


# Convenience function for quick verbose logging
def log_step_verbose(step_id: str, description: Optional[str] = None):
    """Context manager for verbose step logging.
    
    Usage:
        with log_step_verbose("step_033") as logger:
            result = logger.log_operation("Load data", load_data_func)
    """
    class VerboseContext:
        def __init__(self, step_id, description):
            self.logger = VerboseStepLogger(step_id, verbose=True)
            self.description = description
        
        def __enter__(self):
            self.logger.start_step(self.description)
            return self.logger
        
        def __exit__(self, exc_type, exc_val, exc_tb):
            if exc_type:
                self.logger.end_step("error", str(exc_val))
            else:
                self.logger.end_step("success")
    
    return VerboseContext(step_id, description)
