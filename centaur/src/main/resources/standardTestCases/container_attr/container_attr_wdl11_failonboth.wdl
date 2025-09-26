version development-1.1

task dockerAndContainer {
    command <<<
        echo "debian@sha256:9f67f90b1574ea7263a16eb64756897d3fa42a8e43cce61065b8a1f0f9367526"
    >>>
    output {
        File out = stdout()
    }

    # Should fail, these two attributes are mutually exclusive and should not be allowed together
    runtime {
        docker: "ubuntu@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9"
        container: "debian@sha256:9f67f90b1574ea7263a16eb64756897d3fa42a8e43cce61065b8a1f0f9367526"
    }
}

workflow container_attr_wdl11 {
    call dockerAndContainer

    output {
        String out1 = read_string(dockerAndContainer.out)
    }
}
