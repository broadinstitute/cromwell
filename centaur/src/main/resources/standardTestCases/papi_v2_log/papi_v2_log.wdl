version 1.0

workflow papi_v2_log {
    call write_file
    call read_file {
        input: input_file = write_file.written
    }
    call check_log as check_write_log {
        input:
            out_file_path = write_file.out,
            log_file_name = "write_file.log"
    }
    call check_log as check_read_log {
        input:
            out_file_path = read_file.out,
            log_file_name = "read_file.log"
    }
}

task write_file {
    command { echo "That'll do" > written.txt }
    output {
        File out = stdout()
        File written = "written.txt"
    }
    runtime { docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950" }
}

task read_file {
    input { File input_file }
    command { wc -c ~{input_file} }
    output { File out = stdout() }
    runtime { docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950" }
}

task check_log {
    input {
        String out_file_path
        String log_file_name
    }
    String file_log = sub(out_file_path, "/stdout$", "/" + log_file_name)
    command {
        set -euo pipefail
        gsutil cp ~{file_log} log.txt
        set +e
        grep -c 'Starting container setup.' log.txt
        grep -c 'Done container setup.' log.txt
        grep -c 'Starting localization.' log.txt
        grep -c 'Localizing input gs://cloud-cromwell-dev-self-cleaning/cromwell_execution/travis/papi_v2_log/' log.txt
        grep -c 'Done localization.' log.txt
        grep -c 'Running user action: docker run -v /mnt/local-disk:/cromwell_root --entrypoint= ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950 /bin/bash /cromwell_root/script' log.txt
        grep -c 'Starting delocalization.' log.txt
        grep -c 'Delocalizing output /cromwell_root/written.txt -> gs://cloud-cromwell-dev-self-cleaning/cromwell_execution/travis/papi_v2_log/' log.txt
        grep -c 'Delocalizing output /cromwell_root/stdout -> gs://cloud-cromwell-dev-self-cleaning/cromwell_execution/travis/papi_v2_log/' log.txt
        grep -c 'Delocalizing output /cromwell_root/stderr -> gs://cloud-cromwell-dev-self-cleaning/cromwell_execution/travis/papi_v2_log/' log.txt
        grep -c 'Delocalizing output /cromwell_root/rc -> gs://cloud-cromwell-dev-self-cleaning/cromwell_execution/travis/papi_v2_log/' log.txt
        grep -c 'Done delocalization.' log.txt
    }
    output {
        File out = stdout()
        Array[Int] out_counts = read_lines(stdout())
    }
    runtime { docker: "google/cloud-sdk" }
}
