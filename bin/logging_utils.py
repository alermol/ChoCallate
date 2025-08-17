#!/usr/bin/env python3
"""
ChoCallate Python Logging Utilities
Provides structured logging capabilities for Python scripts in the pipeline
"""

import logging
import json
import time
import os
import sys
from datetime import datetime
from typing import Dict, Any, Optional


class ChoCallateLogger:
    """Structured logger for ChoCallate pipeline Python scripts"""
    
    def __init__(self, 
                 script_name: str, 
                 sample_id: Optional[str] = None,
                 log_level: str = "INFO",
                 log_format: str = "json",
                 log_file: Optional[str] = None):
        """
        Initialize the logger
        
        Args:
            script_name: Name of the script
            sample_id: Sample identifier
            log_level: Logging level (DEBUG, INFO, WARN, ERROR, FATAL)
            log_format: Log format (json, text, or both)
            log_file: Optional log file path
        """
        self.script_name = script_name
        self.sample_id = sample_id
        self.log_format = log_format
        self.start_time = time.time()
        
        self.logger = logging.getLogger(script_name)
        self.logger.setLevel(getattr(logging, log_level.upper()))
        
        self.json_formatter = logging.Formatter()
        self.text_formatter = logging.Formatter(
            '%(asctime)s [%(levelname)s] [%(name)s] %(message)s',
            datefmt='%Y-%m-%dT%H:%M:%S'
        )
        
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, log_level.upper()))
        
        if log_format == "json":
            console_handler.setFormatter(self.json_formatter)
        else:
            console_handler.setFormatter(self.text_formatter)
        
        self.logger.addHandler(console_handler)
        
        if log_file:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(getattr(logging, log_level.upper()))
            file_handler.setFormatter(self.json_formatter)
            self.logger.addHandler(file_handler)
        
        self.info("Logger initialized", {
            "script": script_name,
            "sample": sample_id,
            "log_level": log_level,
            "log_format": log_format
        })
    
    def _format_message(self, level: str, message: str, context: Optional[Dict[str, Any]] = None) -> str:
        """Format log message based on configuration"""
        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "level": level,
            "script": self.script_name,
            "sample": self.sample_id,
            "message": message
        }
        
        if context:
            log_entry["context"] = context
        
        if self.log_format == "json":
            return json.dumps(log_entry)
        else:
            parts = [
                f"[{log_entry['timestamp']}]",
                f"[{log_entry['level']}]",
                f"[{log_entry['script']}]"
            ]
            
            if log_entry['sample']:
                parts.append(f"[{log_entry['sample']}]")
            
            parts.append(message)
            
            if context:
                parts.append(f"- Context: {json.dumps(context)}")
            
            return " ".join(parts)
    
    def debug(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Log debug message"""
        formatted_msg = self._format_message("DEBUG", message, context)
        self.logger.debug(formatted_msg)
    
    def info(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Log info message"""
        formatted_msg = self._format_message("INFO", message, context)
        self.logger.info(formatted_msg)
    
    def warn(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Log warning message"""
        formatted_msg = self._format_message("WARN", message, context)
        self.logger.warning(formatted_msg)
    
    def error(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Log error message"""
        formatted_msg = self._format_message("ERROR", message, context)
        self.logger.error(formatted_msg)
    
    def fatal(self, message: str, context: Optional[Dict[str, Any]] = None):
        """Log fatal message"""
        formatted_msg = self._format_message("FATAL", message, context)
        self.logger.critical(formatted_msg)
    
    def log_process_start(self, process_name: str, context: Optional[Dict[str, Any]] = None):
        """Log process start"""
        self.info(f"Process started - {process_name}", context or {})
    
    def log_process_complete(self, process_name: str, context: Optional[Dict[str, Any]] = None):
        """Log process completion"""
        self.info(f"Process completed - {process_name}", context or {})
    
    def log_process_error(self, process_name: str, error: str, context: Optional[Dict[str, Any]] = None):
        """Log process error"""
        error_context = {"error": error}
        if context:
            error_context.update(context)
        self.error(f"Process failed - {process_name}", error_context)
    
    def log_performance(self, metrics: Dict[str, Any]):
        """Log performance metrics"""
        duration = time.time() - self.start_time
        metrics["duration_seconds"] = duration
        self.info("Performance metrics", {"metrics": metrics})
    
    def log_file_operation(self, operation: str, file_path: str):
        """Log file operation"""
        self.info("File operation", {
            "operation": operation,
            "file": file_path
        })
    
    def log_validation(self, validation_type: str, passed: bool, details: Optional[Dict[str, Any]] = None):
        """Log validation result"""
        level = "INFO" if passed else "WARN"
        validation_context = {
            "validation_type": validation_type,
            "passed": passed
        }
        if details:
            validation_context["details"] = details
        
        if level == "INFO":
            self.info(f"Validation {'passed' if passed else 'failed'}", validation_context)
        else:
            self.warn(f"Validation {'passed' if passed else 'failed'}", validation_context)
    
    def get_duration(self) -> float:
        """Get elapsed time since logger initialization"""
        return time.time() - self.start_time


def setup_logger(script_name: str, 
                sample_id: Optional[str] = None,
                log_level: str = "INFO",
                log_format: str = "json",
                log_file: Optional[str] = None) -> ChoCallateLogger:
    """
    Convenience function to set up a logger
    
    Args:
        script_name: Name of the script
        sample_id: Sample identifier
        log_level: Logging level
        log_format: Log format
        log_file: Optional log file path
    
    Returns:
        Configured logger instance
    """
    return ChoCallateLogger(
        script_name=script_name,
        sample_id=sample_id,
        log_level=log_level,
        log_format=log_format,
        log_file=log_file
    )

