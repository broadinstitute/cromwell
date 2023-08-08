task make_file {
    Boolean ready
    command {
        # This comment adds noise to this command to stop it from call caching to other test cases
        # invalidate_bad_caches_use_good
        echo woohoo > out.txt
    }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    }
    output {
        Boolean done = true
        File out = "out.txt"
    }
}

task delete_file_in_gcs {
    Boolean ready
    String file_path
    command {
        gsutil rm ${file_path}
    }
    runtime {
        # google/cloud-sdk:354.0.0-slim
        docker: "google/cloud-sdk@sha256:b5bd0d4b9e56a8b82cea893e7c45f9dfb01fa7cb4e1ce0d426a4468d64654710"
    }
    output {
        Boolean done = true
    }
}

task delete_file_local {
    Boolean ready
    String file_path_raw
    String file_path = sub(file_path_raw, "file://", "")

    command {
        rm ${file_path}
    }
    output {
        Boolean done = true
    }
    runtime {
        backend: "Local"
    }
}

workflow invalidate_bad_caches {
    Boolean running_on_jes

    call make_file as make_first_file { input: ready = true }
    call make_file as make_second_file { input: ready = make_first_file.done }

    if (running_on_jes) {
        call delete_file_in_gcs {
            input:
                ready = make_second_file.done,
                file_path = make_first_file.out
        }
    }
    if (!running_on_jes) {
        call delete_file_local {
            input:
                ready = make_second_file.done,
                file_path_raw = make_first_file.out
        }
    }

    call make_file as cache_third_file {
        input: ready = select_first([delete_file_in_gcs.done, delete_file_local.done])
    }
}
