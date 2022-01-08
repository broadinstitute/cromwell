task make_file {
    command {
        echo woohoo > out.txt
    }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    }
    output {
        File out = "out.txt"
    }
}

task read_file {
    File input_file
    command {
      cat ${input_file}
    }
    runtime {
        docker: "ubuntu:latest"
    }
    output {
        File out = stdout()
    }
}

task delete_file_in_gcs {
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

    call make_file

    if (running_on_jes) {
        call delete_file_in_gcs { input: file_path = make_file.out }
    }
    if (!running_on_jes) {
        call delete_file_local { input: file_path_raw = make_file.out }
    }

    Boolean done = select_first([delete_file_in_gcs.done, delete_file_local.done])

    if (done) {
        call make_file as invalidate_cache_and_remake_file
        call read_file { input: input_file = invalidate_cache_and_remake_file.out }
    }
}

