version development-1.1

workflow InvalidReturnCodeAndContinueOnReturnCode {
    call InvalidReturnCodeContinueOnReturnCode
}

task InvalidReturnCodeContinueOnReturnCode {
    meta {
        volatile: true
    }

    command <<<
        exit 1
    >>>

    runtime {
        docker: "ubuntu:latest"
        returnCodes: [5, 10, 15]
        continueOnReturnCode: [1]
    }
}
