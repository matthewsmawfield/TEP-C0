#!/usr/bin/env python3
"""
TEP-C0 Logging Utilities
================================

Standardized logging infrastructure for the TEP-C0 cosmological analysis pipeline.
Provides color-coded console output, file logging, and custom log levels for
consistent status reporting across all analysis steps.

Custom Log Levels:
    PROCESS (25): Blue - ongoing operations
    SUCCESS (26): Green - completed operations
    TEST (27): Magenta - test/validation results

Author: Matthew Lukin Smawfield
License: CC-BY-4.0
"""

import logging
import sys
from pathlib import Path
import os
from typing import Optional

# Project root for relative path calculations
PACKAGE_ROOT = Path(__file__).resolve().parents[2]

class TEPFormatter(logging.Formatter):
    """Formatter with color support."""

    COLORS = {
        'SUCCESS': '\033[1;32m',  # Green bold
        'WARNING': '\033[1;33m',  # Yellow bold
        'ERROR': '\033[1;31m',    # Red bold
        'INFO': '\033[0;37m',     # White
        'DEBUG': '\033[0;90m',    # Dark gray
        'PROCESS': '\033[0;34m',  # Blue
        'TEST': '\033[1;35m',      # Magenta bold
        'TITLE': '\033[1;36m',     # Cyan bold for titles
        'CRITICAL': '\033[1;41m'  # White on red background
    }
    RESET = '\033[0m'

    def __init__(self, fmt=None, datefmt=None, use_colors=True):
        super().__init__(fmt, datefmt='%H:%M:%S')
        self.use_colors = use_colors

    def format(self, record):
        message = record.getMessage()
        level_mapping = {
            25: ('PROCESS', self.COLORS['PROCESS']),
            26: ('SUCCESS', self.COLORS['SUCCESS']),
            27: ('TEST', self.COLORS['TEST']),
            logging.INFO: ('INFO', self.COLORS['INFO']),
            logging.WARNING: ('WARNING', self.COLORS['WARNING']),
            logging.ERROR: ('ERROR', self.COLORS['ERROR']),
            logging.DEBUG: ('DEBUG', self.COLORS['DEBUG']),
            logging.CRITICAL: ('CRITICAL', self.COLORS['CRITICAL'])
        }

        level_name, color = level_mapping.get(record.levelno, ('INFO', self.COLORS['INFO']))
        timestamp = self.formatTime(record, self.datefmt)
        
        if self.use_colors:
            return f"{color}[{timestamp}] [{level_name}] {message}{self.RESET}"
        else:
            return f"[{timestamp}] [{level_name}] {message}"

class TEPFileFormatter(logging.Formatter):
    """Clean formatter for file output without ANSI color codes."""

    def __init__(self, fmt=None, datefmt=None):
        super().__init__(fmt, datefmt='%H:%M:%S')

    def format(self, record):
        message = record.getMessage()
        level_mapping = {
            25: 'PROCESS',
            26: 'SUCCESS',
            27: 'TEST',
            logging.INFO: 'INFO',
            logging.WARNING: 'WARNING',
            logging.ERROR: 'ERROR',
            logging.DEBUG: 'DEBUG',
            logging.CRITICAL: 'CRITICAL'
        }
        level_name = level_mapping.get(record.levelno, 'INFO')
        timestamp = self.formatTime(record, self.datefmt)
        return f"[{timestamp}] [{level_name}] {message}"

class TEPLogger:
    def __init__(self, name: str = "tep_c0", level: str = "INFO", log_file_path: Optional[Path] = None, reset_log: bool = True):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(self._get_log_level(level))
        self.logger.handlers.clear()
        
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        ch.stream.reconfigure(line_buffering=True)
        ch.setFormatter(TEPFormatter())
        self.logger.addHandler(ch)
        self.logger.propagate = False

        if log_file_path is None:
            default_log_dir = PACKAGE_ROOT / "logs"
            default_log_dir.mkdir(parents=True, exist_ok=True)
            log_file_path = default_log_dir / "tep_c0_pipeline.log"
            reset_log = False

        log_file_path.parent.mkdir(parents=True, exist_ok=True)
        if reset_log:
            try:
                with open(log_file_path, 'w') as f:
                    f.write("")
            except Exception:
                pass
        
        fh = logging.FileHandler(log_file_path, mode='a', encoding='utf-8')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(TEPFileFormatter())
        self.logger.addHandler(fh)

    def _get_log_level(self, level_name: str):
        return getattr(logging, level_name.upper(), logging.INFO)

    def info(self, message: str): self.logger.info(message)
    def warning(self, message: str): self.logger.warning(message)
    def error(self, message: str): self.logger.error(message)
    def debug(self, message: str): self.logger.debug(message)
    def process(self, message: str): self.logger.log(25, message)
    def success(self, message: str): self.logger.log(26, message)
    def test(self, message: str): self.logger.log(27, message)
    def critical(self, message: str): self.logger.critical(message)

_current_step_logger = None

def set_step_logger(logger: TEPLogger):
    global _current_step_logger
    _current_step_logger = logger

def print_status(message: str, level: str = "INFO") -> None:
    if _current_step_logger is not None:
        if level == "SUCCESS": _current_step_logger.success(message)
        elif level == "ERROR": _current_step_logger.error(message)
        elif level == "WARNING": _current_step_logger.warning(message)
        elif level == "PROCESS": _current_step_logger.process(message)
        elif level == "DEBUG": _current_step_logger.debug(message)
        elif level == "TEST": _current_step_logger.test(message)
        elif level == "TITLE":
            _current_step_logger.info("")
            _current_step_logger.info("="*80)
            _current_step_logger.info(f"   {message}")
            _current_step_logger.info("="*80)
            _current_step_logger.info("")
        else: _current_step_logger.info(message)
    else:
        import datetime
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        formatter = TEPFormatter()
        color = formatter.COLORS.get(level, '')
        reset = formatter.RESET
        if level == "TITLE":
            print(f"\n{color}{'='*80}{reset}")
            print(f"{color}   {message}{reset}")
            print(f"{color}{'='*80}{reset}\n")
        else:
            print(f"{color}[{timestamp}] [{level}] {message}{reset}")
    sys.stdout.flush()
