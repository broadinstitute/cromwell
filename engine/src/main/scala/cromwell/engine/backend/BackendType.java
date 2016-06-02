package cromwell.engine.backend;

public enum BackendType {
    LOCAL {
        @Override
        public String displayName() {
            return "Local";
        }
    },
    JES,
    LSF,
    SGE;

    public String displayName() {
        return name();
    }
}
