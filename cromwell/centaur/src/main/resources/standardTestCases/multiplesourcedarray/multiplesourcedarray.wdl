task source {
    Int delay
    command {
        sleep $delay
        echo blah
    }
    output {
        String result = read_string(stdout())
        File resultFile = stdout()
    }
    runtime {docker: "ubuntu:latest"}
}

task consume {
    Array[String] in_file
    command {
        echo ${sep=' ' in_file}
    }
    output {
        String out = read_string(stdout())
    }
    runtime {docker: "ubuntu:latest"}
}

task consumeFiles {
    Array[File] in_file
    command {
        cat ${sep='; cat ' in_file}
    }
    output {
        String x = read_string(stdout())
    }
    runtime {docker: "ubuntu:latest"}
}

workflow multplesourcedarray {
    call source as source1 { input: delay=1 }
    call source as source2 { input: delay=2 }
    call source as source3 { input: delay=3 }
    call source as source4 { input: delay=4 }
    call source as source5 { input: delay=10 }
    call source as source6 { input: delay=11 }
    call source as source7 { input: delay=12 }

    call consume { input: in_file = [source1.result, source2.result, source3.result, source4.result, source5.result, source6.result, source7.result] }
    call consumeFiles { input: in_file = [source1.resultFile, source2.resultFile, source3.resultFile, source4.resultFile, source5.resultFile, source6.resultFile, source7.resultFile] }

    output {
        consumeFiles.x
    }
}
