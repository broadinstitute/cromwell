package cromwell.parser;

import static cromwell.parser.BackendType.JES;
import static cromwell.parser.BackendType.LOCAL;

/**
 * Backend runtime keys and the backends which are known to support them.
 */
public enum RuntimeKey {
    CPU("cpu", JES),
    DEFAULT_DISKS("defaultDisks", JES),
    DEFAULT_ZONES("defaultZones", JES),
    DOCKER("docker", new BackendType[]{JES}, LOCAL),
    FAIL_ON_STDERR("failOnStderr", LOCAL, JES),
    MEMORY("memory", JES),
    PREEMPTIBLE("preemptible", JES);

    public final String key;

    /**
     * BackendTypes on which this key is mandatory.
     */
    public final BackendType[] mandatory;

    /**
     * BackendTypes optionally supporting this key.  A BackendType need only appear in one
     * of the mandatory or optional arrays.
     */
    public final BackendType[] optional;

    RuntimeKey(String key, BackendType... optional) {
        this(key, new BackendType[0], optional);
    }

    RuntimeKey(String key, BackendType [] mandatory, BackendType... optional) {
        this.key = key;
        this.mandatory = mandatory;
        this.optional = optional;
    }

    public boolean isMandatory(BackendType backendType) {
        for (BackendType type : mandatory) {
            if (type == backendType) {
                return true;
            }
        }
        return false;
    }

    public boolean isOptional(BackendType backendType) {
        for (BackendType type : optional) {
            if (type == backendType) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns true if the specified key is mandatory or optional on this backend.
     */
    public boolean supports(BackendType backendType) {
        return isOptional(backendType) || isMandatory(backendType);
    }
}
