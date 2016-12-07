task output_redirection {
    command {
        echo "#!/bin/bash" > writeToStderr.sh
        # Yes, this is deliberate. We'll eventually end up redirecting this to stdout:
        echo "echo 'should be on stdout' >&2" >> writeToStderr.sh

        chmod +x writeToStderr.sh

        # Write to stdout but redirect to stderr:
        echo "should be on stderr" >&2
        # Write to stderr but redirect to stdout:
        ./writeToStderr.sh 2>&1

        sleep 2
    }
    output {
        String stdout = read_string(stdout())
        String stderr = read_string(stderr())
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

workflow output_redirection_wf {
    call output_redirection
}
