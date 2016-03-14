package cromwell.engine.backend;

public enum BackendType {
    LOCAL {
        @Override
        public String displayName() {
            return "Local";
        }
    },
    JES,
    PBS,
    SGE;

    public String displayName() {
        return name();
    }
}
