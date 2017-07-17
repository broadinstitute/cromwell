package wdl4s.parser;

public enum MemoryUnit {
    Bytes(1, "B"),
    KiB(1L << 10, "KiB", "Ki"),
    MiB(1L << 20, "MiB", "Mi"),
    GiB(1L << 30, "GiB", "Gi"),
    TiB(1L << 40, "TiB", "Ti"),
    KB(1000, "KB", "K"),
    MB(KB.bytes * 1000, "MB", "M"),
    GB(MB.bytes * 1000, "GB", "G"),
    TB(GB.bytes * 1000, "TB", "T");

    public final double bytes;
    public final String[] suffixes;

    MemoryUnit(double bytes, String... suffixes) {
        this.bytes = bytes;
        this.suffixes = suffixes;
    }

    public static MemoryUnit fromSuffix(String suffix) {
        for (MemoryUnit unit: values()) {
            for (String unitSuffix: unit.suffixes) {
                if (unitSuffix.equals(suffix)) return unit;
            }
        }

        throw new IllegalArgumentException("Unit with suffix " + suffix + " was not found.");
    }
}

