package cromwell.parser;

import com.typesafe.config.Config;

public enum BackendType {
    LOCAL, JES, SGE;
    public static BackendType from(Config backendConf) {
        return from(backendConf.getString("backend"));
    }

    public static BackendType from(String name) {
        switch (name.toLowerCase()) {
            case "local":
                return LOCAL;
            case "jes":
                return JES;
            case "sge":
                return SGE;
            default:
                throw new IllegalArgumentException(name + " is not a recognized backend");
        }
    }
}
