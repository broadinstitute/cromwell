task dockerOnly {
    command <<<
        echo "Run with WDL draft-2 on ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: "ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
    }
}

task dockerAndContainer {
    command <<<
        echo "Run with WDL draft-2 on ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: "ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
        container: "debian@sha256:9f67f90b1574ea7263a16eb64756897d3fa42a8e43cce61065b8a1f0f9367526"
    }
}

workflow container_attr_wdldraft2 {
    call dockerOnly
    call dockerAndContainer

    output {
        String out1 = read_string(dockerOnly.out)
        String out2 = read_string(dockerAndContainer.out)
    }
}
