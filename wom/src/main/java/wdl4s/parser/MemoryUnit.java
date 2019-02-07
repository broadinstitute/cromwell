package wdl4s.parser;

public enum MemoryUnit {
    Bytes(1, "B"),
    KB(1L << 10, "KB", "K", "KiB", "Ki"),
    MB(1L << 20, "MB", "M", "MiB", "Mi"),
    GB(1L << 30, "GB", "G", "GiB", "Gi"),
    TB(1L << 40, "TB", "T", "TiB", "Ti");

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

