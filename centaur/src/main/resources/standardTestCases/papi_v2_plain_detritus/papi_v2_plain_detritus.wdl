version 1.0

workflow papi_v2_plain_detritus {
    call get_stdout
    call check_detritus {
        input:
            out_file_path = get_stdout.out
    }
}

task get_stdout {
    command { echo "Hello stdout" }
    output {
        File out = stdout()
    }
    runtime { docker: "google/cloud-sdk:251.0.0-slim" }
}


task check_detritus {
    input {
        String out_file_path
    }
    String rc_file = sub(out_file_path, "/stdout$", "/rc")
    command {
        gsutil ls -L ~{rc_file} | grep 'Content-Type:' | sed 's/.* Content-Type: *//' | sed 's/ *$//'
    }
    output {
        String out = read_string(stdout())
    }
    runtime { docker: "google/cloud-sdk:251.0.0-slim" }
}
