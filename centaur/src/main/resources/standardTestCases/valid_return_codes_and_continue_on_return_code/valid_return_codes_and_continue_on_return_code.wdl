version development-1.1

workflow ValidReturnCodeAndContinueOnReturnCode {
    call ReturnCodeContinueOnReturnCode1
    call ReturnCodeContinueOnReturnCode2
    call ReturnCodeContinueOnReturnCode3
}

task ReturnCodeContinueOnReturnCode1 {
    meta {
        volatile: true
    }

    command <<<
        exit 1
    >>>

    runtime {
        docker: "ubuntu:latest"
        returnCodes: [1]
        continueOnReturnCode: [0]
    }
}

task ReturnCodeContinueOnReturnCode2 {
    meta {
        volatile: true
    }

    command <<<
        exit 1
    >>>

    runtime {
        docker: "ubuntu:latest"
        returnCodes: [1]
        continueOnReturnCode: false
    }
}

task ReturnCodeContinueOnReturnCode3 {
    meta {
        volatile: true
    }

    command <<<
        exit 1
    >>>

    runtime {
        docker: "ubuntu:latest"
        returnCodes: [1, 4, 7]
        continueOnReturnCode: [1, 3, 5]
    }
}
