version development-1.1

task dockerSingle {
    command <<<
        echo "Run on ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: "ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
    }
}

task containerSingle {
    command <<<
        echo "Run on debian@sha256:9f67f90b1574ea7263a16eb64756897d3fa42a8e43cce61065b8a1f0f9367526"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        container: "debian@sha256:9f67f90b1574ea7263a16eb64756897d3fa42a8e43cce61065b8a1f0f9367526"
    }
}

task dockerList {
    command <<<
        echo "Run on ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: [
            "ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9",
            "debian@sha256:9f67f90b1574ea7263a16eb64756897d3fa42a8e43cce61065b8a1f0f9367526"
        ]
    }
}

task containerList {
    command <<<
        echo "Run on debian@sha256:9f67f90b1574ea7263a16eb64756897d3fa42a8e43cce61065b8a1f0f9367526"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        container: [
           "debian@sha256:9f67f90b1574ea7263a16eb64756897d3fa42a8e43cce61065b8a1f0f9367526",
           "ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
        ]
    }
}

workflow container_attr_wdl11 {
    call dockerSingle
    call containerSingle
    call dockerList
    call containerList

    output {
        String out1 = read_string(dockerSingle.out)
        String out2 = read_string(containerSingle.out)
        String out3 = read_string(dockerList.out)
        String out4 = read_string(containerList.out)
    }
}
