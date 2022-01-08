version 1.0

workflow workbench_health_monitor_check {
    call run
    output {
        String out = run.out
    }
}

task run {
    command {
        # Try to apt-get install jq and curl if they are missing. These should already be present in Travis so this would be a noop.
        # Do this all quietly with output redirected to /dev/null since the logic below reads stdout.
        (which jq > /dev/null && which curl > /dev/null) || (apt-get update > /dev/null && apt-get install -y jq curl > /dev/null)

        # If a Cromwell port has been specified on the command line use that, otherwise assume the managed Cromwell is on 8008.
        port=$(ps auxwww | grep /app/cromwell.jar | grep -v grep | sed -E 's/.*-Dwebservice\.port=([0-9]+).*/\1/')
        if [ -z "$port" ]; then port=8008; fi

        set -euo pipefail
        curl http://localhost:$port/engine/v1/status | jq -crMS .
    }
    runtime { backend: "Local" }
    output { String out = read_string(stdout()) }
}
