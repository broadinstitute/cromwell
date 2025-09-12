version 1.0

task dockerOnly {
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

task dockerAndContainer {
    command <<<
        echo "Run on ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: "ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
        container: "alpine:latest"
    }
}

task dockerAndContainerList {
    command <<<
        echo "Run on ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: "ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
        container: ["alpine:latest", "debian:latest"]
    }
}

workflow container_attr_wdl10 {
    call dockerOnly
    call dockerAndContainer
    call dockerAndContainerList

    output {
        String out1 = read_string(dockerOnly.out)
        String out2 = read_string(dockerAndContainer.out)
        String out3 = read_string(dockerAndContainerList.out)
    }
}
