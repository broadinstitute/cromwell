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
    DOCKER("docker", LOCAL, JES),
    FAIL_ON_STDERR("failOnStderr", LOCAL, JES),
    MEMORY("memory", JES),
    PREEMPTIBLE("preemptible", JES);

    public final String key;

    public final BackendType[] backendTypes;

    RuntimeKey(String key, BackendType... backendTypes) {
        this.key = key;
        this.backendTypes = backendTypes;
    }
}
