task make_file {
    Int content
    command {
        echo "hello ${content}" > out
    }
    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:53a002b59dfcd43b4d15e97c1acbeae035ddd1b31a106659a312e9fe65f00afa"
        backend: "Papi-Caching-No-Copy"
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
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:53a002b59dfcd43b4d15e97c1acbeae035ddd1b31a106659a312e9fe65f00afa"
        backend: "Papi-Caching-No-Copy"
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
