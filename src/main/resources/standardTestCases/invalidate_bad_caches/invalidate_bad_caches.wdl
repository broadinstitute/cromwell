task make_file {
    command {
        echo woohoo > out.txt && sleep 2
    }
    runtime {
        docker: "ubuntu:latest"
    }
    output {
        File out = "out.txt"
    }
}

task delete_file_in_gcs {
    String file_path
    command {
        gsutil rm ${file_path} && sleep 2
    }
    runtime {
        docker: "google/cloud-sdk"
    }
    output {
        Boolean done = true
    }
}

task delete_file_local {
    String file_path_raw
    String file_path = sub(file_path_raw, "file://", "")

    command {
        rm ${file_path} && sleep 2
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
    }
}

