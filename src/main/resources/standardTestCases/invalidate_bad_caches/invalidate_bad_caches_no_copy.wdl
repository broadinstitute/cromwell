task make_file {
    Boolean ready = true
    command {
        echo woohoo > out.txt
    }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        backend: "Jes-Caching-No-Copy"
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
        backend: "Jes-Caching-No-Copy"
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
        docker: "google/cloud-sdk@sha256:fb904276e8a902ccd9564989d9222bdfbe37ffcd7f9989ca7e24b4019a9b4b6b"
        backend: "Jes-Caching-No-Copy"
    }
    output {
        Boolean done = true
    }
}

workflow invalidate_bad_caches {
    call make_file

    call delete_file_in_gcs { input: file_path = make_file.out }

    call make_file as invalidate_cache_and_remake_file { input: ready = delete_file_in_gcs.done }
    call read_file { input: input_file = invalidate_cache_and_remake_file.out }
}
