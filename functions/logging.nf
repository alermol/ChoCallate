// ChoCallate Structured Logging Module
// Provides consistent logging format and levels across all pipeline processes

import groovy.json.JsonBuilder
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

// Log levels
enum LogLevel {
    DEBUG(0), INFO(1), WARN(2), ERROR(3), FATAL(4)
    
    private final int level
    
    LogLevel(int level) {
        this.level = level
    }
    
    int getLevel() { return level }
}

// Get current timestamp in ISO format
def getCurrentTimestamp() {
    return LocalDateTime.now().format(DateTimeFormatter.ISO_LOCAL_DATE_TIME)
}

// Get logging configuration - creates a new config object each time
def getLogConfig() {
    try {
        def outdir = (params && params.outdir) ? params.outdir : 'ChoCallate_output'
        def logFile = (params && params.log_file) ? params.log_file : 'ChoCallate.log'
        def errorLogFile = (params && params.log_error_file) ? params.log_error_file : 'ChoCallate_errors.log'
        
        // Construct full paths if outdir is specified
        if (outdir && outdir != '.') {
            logFile = "${outdir}/${logFile}"
            errorLogFile = "${outdir}/${errorLogFile}"
        }
        
        return [
            level: (params && params.log_level) ? params.log_level : 'INFO',
            format: (params && params.log_format) ? params.log_format : 'json',
            include_timestamp: (params && params.log_timestamp != null) ? params.log_timestamp : true,
            include_process: (params && params.log_process != null) ? params.log_process : true,
            include_sample: (params && params.log_sample != null) ? params.log_sample : true,
            log_file: logFile,
            error_log_file: errorLogFile
        ]
    } catch (Exception e) {
        // Fallback to defaults if params access fails
        return [
            level: 'INFO',
            format: 'json',
            include_timestamp: true,
            include_process: true,
            include_sample: true,
            log_file: 'ChoCallate.log',
            error_log_file: 'ChoCallate_errors.log'
        ]
    }
}

// Check if logging level should be displayed
def shouldLog(LogLevel level) {
    def config = getLogConfig()
    def configLevel = LogLevel.valueOf(config.level.toUpperCase())
    return level.getLevel() >= configLevel.getLevel()
}

// Core logging function
def logMessage(LogLevel level, String message, Map context = [:], String processName = null, String sampleId = null) {
    if (!shouldLog(level)) return
    
    def logEntry = [
        timestamp: getCurrentTimestamp(),
        level: level.name(),
        message: message,
        context: context
    ]
    
    def config = getLogConfig()
    if (config.include_process && processName) {
        logEntry.process = processName
    }
    
    if (config.include_sample && sampleId) {
        logEntry.sample = sampleId
    }
    
    // Add workflow context if available
    try {
        if (workflow) {
            logEntry.workflow = [
                runName: workflow.runName,
                sessionId: workflow.sessionId,
                projectDir: workflow.projectDir?.toString()
            ]
        }
    } catch (Exception e) {
        // Workflow context not available
    }
    
    // Format and output log entry
    def formattedLog = formatLogEntry(logEntry, config)
    
    // Write to appropriate log file
    writeToLogFile(logEntry, level, config)
    
    // Also output to console for Nextflow processes
    if (level == LogLevel.ERROR || level == LogLevel.FATAL) {
        println "ERROR: ${message}"
        if (context) {
            println "Context: ${new JsonBuilder(context).toString()}"
        }
    } else if (level == LogLevel.WARN) {
        println "WARNING: ${message}"
    } else if (level == LogLevel.INFO) {
        println "INFO: ${message}"
    } else if (level == LogLevel.DEBUG) {
        try {
            if (params && params.debug) {
                println "DEBUG: ${message}"
            }
        } catch (Exception e) {
            // Debug mode not available
        }
    }
}

// Format log entry based on configuration
def formatLogEntry(Map logEntry, Map config) {
    if (config.format == 'json') {
        return new JsonBuilder(logEntry).toString()
    } else if (config.format == 'text') {
        return formatTextLog(logEntry)
    } else {
        // Both formats
        return "${formatTextLog(logEntry)}\n${new JsonBuilder(logEntry).toString()}"
    }
}

// Format text log entry
def formatTextLog(Map logEntry) {
    def parts = []
    parts << "[${logEntry.timestamp}]"
    parts << "[${logEntry.level}]"
    
    if (logEntry.process) {
        parts << "[${logEntry.process}]"
    }
    
    if (logEntry.sample) {
        parts << "[${logEntry.sample}]"
    }
    
    parts << logEntry.message
    
    if (logEntry.context) {
        parts << "- Context: ${logEntry.context.toString()}"
    }
    
    return parts.join(" ")
}

// Write log entry to appropriate file
def writeToLogFile(Map logEntry, LogLevel level, Map config) {
    try {
        def logFile = (level == LogLevel.ERROR || level == LogLevel.FATAL) ? 
                     config.error_log_file : config.log_file
        
        def file = new File(logFile)
        if (!file.exists()) {
            def parentDir = file.getParentFile()
            if (parentDir && !parentDir.exists()) {
                parentDir.mkdirs()
            }
            file.createNewFile()
        }
        
        file.append("${new JsonBuilder(logEntry).toString()}\n")
    } catch (Exception e) {
        // Fallback to console if file writing fails
        println "Failed to write to log file: ${e.message}"
    }
}

// Convenience logging functions
def logDebug(String message, Map context = [:], String processName = null, String sampleId = null) {
    logMessage(LogLevel.DEBUG, message, context, processName, sampleId)
}

def logInfo(String message, Map context = [:], String processName = null, String sampleId = null) {
    logMessage(LogLevel.INFO, message, context, processName, sampleId)
}

def logWarn(String message, Map context = [:], String processName = null, String sampleId = null) {
    logMessage(LogLevel.WARN, message, context, processName, sampleId)
}

def logError(String message, Map context = [:], String processName = null, String sampleId = null) {
    logMessage(LogLevel.ERROR, message, context, processName, sampleId)
}

def logFatal(String message, Map context = [:], String processName = null, String sampleId = null) {
    logMessage(LogLevel.FATAL, message, context, processName, sampleId)
}

// Process-specific logging functions
def logProcessStart(String processName, String sampleId, Map context = [:]) {
    logInfo("Process started", context + [action: "start"], processName, sampleId)
}

def logProcessComplete(String processName, String sampleId, Map context = [:]) {
    logInfo("Process completed", context + [action: "complete"], processName, sampleId)
}

def logProcessError(String processName, String sampleId, String error, Map context = [:]) {
    logError("Process failed", context + [action: "error", error: error], processName, sampleId)
}

// Workflow lifecycle logging
def logWorkflowStart() {
    try {
        if (workflow) {
            logInfo("Workflow execution started", [
                action: "workflow_start",
                runName: workflow.runName,
                sessionId: workflow.sessionId
            ])
        } else {
            logInfo("Workflow execution started", [
                action: "workflow_start",
                note: "Workflow context not available"
            ])
        }
    } catch (Exception e) {
        logInfo("Workflow execution started", [
            action: "workflow_start",
            note: "Workflow context not available"
        ])
    }
}

def logWorkflowComplete() {
    try {
        if (workflow) {
            logInfo("Workflow execution completed", [
                action: "workflow_complete",
                runName: workflow.runName,
                sessionId: workflow.sessionId
            ])
        } else {
            logInfo("Workflow execution completed", [
                action: "workflow_complete",
                note: "Workflow context not available"
            ])
        }
    } catch (Exception e) {
        logInfo("Workflow execution completed", [
            action: "workflow_complete",
            note: "Workflow context not available"
        ])
    }
}

def logWorkflowError(String error) {
    try {
        if (workflow) {
            logFatal("Workflow execution failed", [
                action: "workflow_error",
                error: error,
                runName: workflow.runName,
                sessionId: workflow.sessionId
            ])
        } else {
            logFatal("Workflow execution failed", [
                action: "workflow_error",
                error: error,
                note: "Workflow context not available"
            ])
        }
    } catch (Exception e) {
        logFatal("Workflow execution failed", [
            action: "workflow_error",
            error: error,
            note: "Workflow context not available"
        ])
    }
}

// Performance logging
def logPerformance(String processName, String sampleId, Map metrics) {
    logInfo("Performance metrics", [
        action: "performance",
        metrics: metrics
    ], processName, sampleId)
}

// Resource usage logging
def logResourceUsage(String processName, String sampleId, Map resources) {
    logInfo("Resource usage", [
        action: "resource_usage",
        resources: resources
    ], processName, sampleId)
}

// File operation logging
def logFileOperation(String operation, String filePath, String processName = null, String sampleId = null) {
    logInfo("File operation", [
        action: "file_operation",
        operation: operation,
        file: filePath
    ], processName, sampleId)
}

// Validation logging
def logValidation(String validationType, boolean passed, Map details = [:], String processName = null, String sampleId = null) {
    def level = passed ? LogLevel.INFO : LogLevel.WARN
    logMessage(level, "Validation ${passed ? 'passed' : 'failed'}", [
        action: "validation",
        type: "validation",
        passed: passed,
        details: details
    ], processName, sampleId)
}

// Initialize logging - creates log directory and logs pipeline start
def initLogging() {
    def config = getLogConfig()
    
    // Create log directory if it doesn't exist
    def logDir = new File(config.log_file).getParentFile()
    if (logDir && !logDir.exists()) {
        logDir.mkdirs()
    }
    // If no parent directory (relative path), logs will be created in current directory
    
    // Log pipeline start - only if params are available
    try {
        if (params) {
            logInfo("Pipeline started", [
                pipeline: "ChoCallate",
                version: "1.0.0",
                timestamp: getCurrentTimestamp(),
                parameters: [
                    samples_tsv: params.samples_tsv,
                    outdir: params.outdir,
                    ploidy: params.ploidy,
                    reads_type: params.reads_type,
                    reads_source: params.reads_source,
                    effective_callers: params.effective_callers,
                    min_coverage: params.min_coverage,
                    min_base_quality: params.min_base_quality,
                    min_map_qual: params.min_map_qual,
                    min_snp_qual: params.min_snp_qual,
                    cons_type: params.cons_type,
                    debug: params.debug
                ]
            ])
        } else {
            logInfo("Pipeline started", [
                pipeline: "ChoCallate",
                version: "1.0.0",
                timestamp: getCurrentTimestamp(),
                note: "Parameters not yet available"
            ])
        }
    } catch (Exception e) {
        logInfo("Pipeline started", [
            pipeline: "ChoCallate",
            version: "1.0.0",
            timestamp: getCurrentTimestamp(),
            note: "Parameters not available due to error: ${e.message}"
        ])
    }
}
