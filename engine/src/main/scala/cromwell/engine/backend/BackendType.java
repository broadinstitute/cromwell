package cromwell.engine.backend;

@Deprecated
public enum BackendType {
    LOCAL {
        @Override
        public String displayName() {
            return "Local";
        }
    },
    JES,
    SGE;

    public String displayName() {
        return name();
    }
}
