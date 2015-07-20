package cromwell.parser;

/**
 * Convenience enum for parsing memory specifications.  If feels like there should be a preexisting library
 * for this but I couldn't find one.
 */
public enum MemorySize {
    Bytes(1, "B"),
    KiB(1_000, "KiB"),
    MiB(KiB.multiplier * KiB.multiplier, "MiB"),
    GiB(KiB.multiplier * MiB.multiplier, "GiB"),
    TiB(KiB.multiplier * GiB.multiplier, "TiB"),
    KB(1 << 10, "KB"),
    MB(KB.multiplier * KB.multiplier, "MB"),
    GB(KB.multiplier * MB.multiplier, "GB"),
    TB(KB.multiplier * GB.multiplier, "TB");

    public final long multiplier;
    public final String suffix;

    MemorySize(long multiplier, String suffix) {
        this.multiplier = multiplier;
        this.suffix = suffix;
    }

    /** Convert from the units of this memory size to individual bytes. */
    public double toBytes(double prefix) {
        return prefix * multiplier;
    }

    /** Convert from individual bytes to units of this memory size. */
    public double fromBytes(double bytes) {
        return bytes / multiplier;
    }
}
