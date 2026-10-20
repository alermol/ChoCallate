/*
 * Helpers shared by process scripts to size JVM heaps from task resources.
 */
class ResourceUtils {

    /**
     * Heap size in MiB for one JVM.
     *
     * A task may start several JVMs in parallel (GNU parallel), in which
     * case each receives an equal share.  Twenty percent of the task
     * memory is left for native allocations, stacks and the file cache.
     *
     * @param memory task memory, as exposed by Nextflow
     * @param slots number of JVMs that may run concurrently
     * @return heap size in MiB, or 0 when no memory is configured
     */
    static int javaHeapMb(memory, int slots = 1) {
        if (memory == null) {
            return 0
        }
        int heapSlots = slots > 0 ? slots : 1
        return (memory.toMega() * 0.8 / heapSlots) as int
    }

    /**
     * GATK/Picard {@code -Xmx} option for one JVM.
     */
    static String javaHeapOption(memory, int slots = 1) {
        int heapMb = javaHeapMb(memory, slots)
        return heapMb > 0 ? "-Xmx${heapMb}m" : ''
    }
}
