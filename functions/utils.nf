// Function to get the number of available callers
def getAvailableCallersCount(String available_callers) {
    def callersList = available_callers.split(/\s*,\s*/).findAll { it }
    return callersList.size()
}

// Function to get the consensus threshold based on the consensus type
def getConsensusThreshold(String cons_type, String available_callers) {
    if (cons_type == 'mj') {
        return getAvailableCallersCount(available_callers).intdiv(2) + 1
    } else if (cons_type == 'n1') {
        return getAvailableCallersCount(available_callers) - 1
    } else if (cons_type == 'fc') {
        return getAvailableCallersCount(available_callers)
    } else {
        error "Invalid consensus type: ${cons_type}. Available types: mj, n1, fc"
    }
}