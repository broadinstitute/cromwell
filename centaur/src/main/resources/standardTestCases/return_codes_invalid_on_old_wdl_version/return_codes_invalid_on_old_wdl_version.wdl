workflow ReturnCodesInvalidOnOldWdl {
        call ReturnCodesInvalidOnOldWdlTask
}

task ReturnCodesInvalidOnOldWdlTask {
    meta {
        volatile: "true"
    }

    command <<<
        exit 5
    >>>

    runtime {
        docker: "ubuntu:latest"
        returnCodes: [5, 10, 15]
    }
}
