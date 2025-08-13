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

// Function to check if all effective callers are in the available callers
def allEffectiveCallersInAvailable(String effective_callers, String available_callers) {
    def effList = effective_callers.split(/\s*,\s*/).findAll { it }
    def availList = available_callers.split(/\s*,\s*/).findAll { it }
    def missing = effList.findAll { !(it in availList) }
    if (missing) {
        println "The following effective callers are not in available callers: ${missing.join(', ')}"
        return false
    }
    return true
}

// Function to check if the number of effective callers is at least three
def effectiveCallersAtLeastThree(String effective_callers) {
    def effList = effective_callers.split(/\s*,\s*/).findAll { it }
    return effList.size()
}

// Function to check if all effective callers are suitable for diploid calling when ploidy == 2
def allEffectiveCallersDiploidSuitable(String effective_callers, String diploid_callers) {
    def effList = effective_callers.split(/\s*,\s*/).findAll { it }
    def dipList = diploid_callers.split(/\s*,\s*/).findAll { it }
    def unsuitable = effList.findAll { !(it in dipList) }
    if (unsuitable) {
        println "The following effective callers are not suitable for diploid calling: ${unsuitable.join(', ')}"
        return false
    }
    return true
}

// Function to check if all effective callers are suitable for polyploid calling
def allEffectiveCallersPolyploidSuitable(String effective_callers, String polyploid_callers) {
    def effList = effective_callers.split(/\s*,\s*/).findAll { it }
    def polyList = polyploid_callers.split(/\s*,\s*/).findAll { it }
    def unsuitable = effList.findAll { !(it in polyList) }
    if (unsuitable) {
        println "The following effective callers are not suitable for polyploid calling: ${unsuitable.join(', ')}"
        return false
    }
    return true
}

