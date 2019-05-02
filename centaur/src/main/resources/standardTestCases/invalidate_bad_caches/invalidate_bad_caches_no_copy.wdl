task make_file {
    Boolean ready = true
    command {
        # This comment adds noise to this command to stop it from call caching to other test cases
        echo woohoo > out.txt
    }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        backend: "Papi-Caching-No-Copy"
    }
    output {
        Boolean done = true
        File out = "out.txt"
    }
}

task read_file {
    File input_file
    command {
      # This comment adds noise to this command to stop it from call caching to other test cases
      cat ${input_file}
    }
    runtime {
        docker: "ubuntu:latest"
        backend: "Papi-Caching-No-Copy"
    }
    output {
        String out = read_string(stdout())
    }
}

task delete_file_in_gcs {
    String file_path
    command {
        # This comment adds noise to this command to stop it from call caching to other test cases
        gsutil rm ${file_path}
    }
    runtime {
        docker: "google/cloud-sdk@sha256:fb904276e8a902ccd9564989d9222bdfbe37ffcd7f9989ca7e24b4019a9b4b6b"
        backend: "Papi-Caching-No-Copy"
    }
    output {
        Boolean done = true
    }
}

workflow invalidate_bad_caches_no_copy {
    # Make a file the first time:
    call make_file

    # Because it will call cache, we'll get the same file here.
    call make_file as make_file_again { input: ready = make_file.done }

    # Delete both because we referenced (rather than copied) the file
    call delete_file_in_gcs { input: file_path = make_file_again.out }

    # Re-make the file:
    call make_file as invalidate_cache_and_remake_file { input: ready = delete_file_in_gcs.done }

    call read_file { input: input_file = invalidate_cache_and_remake_file.out }

    output {
      String woohoo = read_file.out
    }
}
