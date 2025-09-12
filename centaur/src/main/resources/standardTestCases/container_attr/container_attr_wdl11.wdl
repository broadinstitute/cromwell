version 1.1

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

task containerOnly {
    command <<<
        echo "Run on alpine@sha256:4bcff63911fcb4448bd4fdacec207030997caf25e9bea4045fa6c8c44de311d1"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        container: "alpine@sha256:4bcff63911fcb4448bd4fdacec207030997caf25e9bea4045fa6c8c44de311d1"
    }
}

task dockerAndContainer {
    command <<<
        echo "Run on alpine@sha256:4bcff63911fcb4448bd4fdacec207030997caf25e9bea4045fa6c8c44de311d1"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: "ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
        container: "alpine@sha256:4bcff63911fcb4448bd4fdacec207030997caf25e9bea4045fa6c8c44de311d1"
    }
}

task dockerAndContainerList {
    command <<<
        echo "Run on alpine@sha256:4bcff63911fcb4448bd4fdacec207030997caf25e9bea4045fa6c8c44de311d1"
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: "ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
        container: [
            "alpine@sha256:4bcff63911fcb4448bd4fdacec207030997caf25e9bea4045fa6c8c44de311d1",
            "debian@sha256:9f67f90b1574ea7263a16eb64756897d3fa42a8e43cce61065b8a1f0f9367526"
        ]
    }
}

workflow runtime_continueOnRC {
    call dockerOnly
    call containerOnly
    call dockerAndContainer
    call dockerAndContainerList

    output {
        String out1 = read_string(dockerOnly.out)
        String out2 = read_string(containerOnly.out)
        String out3 = read_string(dockerAndContainer.out)
        String out4 = read_string(dockerAndContainerList.out)
    }
}
