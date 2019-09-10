workflow Largenumber {
    Float largenumber

    call print_number {
        input: largenumber=largenumber
    }
}
task print_number {
    Float largenumber
    Float largeNumIncr = largenumber + 42
    command {
        echo SOME_LARGE_NUMBER  ${largeNumIncr}
    }

    runtime {
      docker: "ubuntu"
      memory: "2G"
      cpu: 1
    }

    output {
        String lnum=read_string(stdout())
    }
}
