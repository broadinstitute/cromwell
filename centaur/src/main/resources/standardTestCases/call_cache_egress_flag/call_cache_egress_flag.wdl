workflow call_cache_egress_flag {
    # Make some files:
    call make_files { input: ready = true }

    # Call cache the made files
    call make_files as make_files_cached { input: ready = make_files.done }
}

task make_files {
    Boolean ready
    command {
        echo bananeira > bananeira.txt
        echo balanca > balanca.txt
    }
    runtime { docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950" }
    output {
        Boolean done = true
        File bananeira = "bananeira.txt"
        File balanca = "balanca.txt"
    }
}