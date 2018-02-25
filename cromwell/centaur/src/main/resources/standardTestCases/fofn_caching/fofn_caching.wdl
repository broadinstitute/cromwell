task make_file {
    Int content
    command {
        echo "hello ${content}" > out
    }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        backend: "Jes-Caching-No-Copy"
    }
    output {
        File out = "out"
    }
}

task use_fofn {
    File fofn
    command {
        cat ${fofn}
    }
    runtime {
        docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
        backend: "Jes-Caching-No-Copy"
    }
    output {
        Array[String] out = read_lines(stdout())
    }
}


workflow fofn_caching {
    scatter(i in range(5)) {
        call make_file { input: content = i }
    }
    
    call use_fofn { input: fofn = write_lines(make_file.out) }
    
    output {
        Array[String] files = use_fofn.out
    }
}
