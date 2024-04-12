version development-1.1

workflow FailingReturnCode {
    call FailingReturnCodeSet
}

task FailingReturnCodeSet {
    meta {
        volatile: true
    }

    command <<<
        exit 1
    >>>

    runtime {
        returnCodes: [0, 5, 10]
        docker: "ubuntu:latest"
    }
}
