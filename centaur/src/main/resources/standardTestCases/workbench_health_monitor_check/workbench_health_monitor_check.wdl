version 1.0

workflow workbench_health_monitor_check {
    call run
    output {
        String out = run.out
    }
}

task run {
    command {
        set -euo pipefail
        curl http://localhost:8008/engine/v1/status | jq -crMS .
    }
    runtime { backend: "Local" }
    output { String out = read_string(stdout()) }
}
