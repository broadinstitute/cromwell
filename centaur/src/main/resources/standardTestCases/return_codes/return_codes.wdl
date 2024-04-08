version development-1.1

workflow ReturnCodeValidation {
    call ReturnCodeSet1
    call ReturnCodeSet2
    call ReturnCodeSet3
    call ReturnCodeString
}

task ReturnCodeSet1 {
    meta {
        volatile: true
    }

    command <<<
        exit 1
    >>>

    runtime {
        docker: "ubuntu:latest"
        returnCodes: [1]
    }
}

task ReturnCodeSet2 {
    meta {
        volatile: true
    }

    command <<<
        exit 200
    >>>

    runtime {
        docker: "ubuntu:latest"
        returnCodes: [1, 123, 200]
    }
}

task ReturnCodeSet3 {
    meta {
        volatile: true
    }

    command <<<
        exit 10
    >>>

    runtime {
        docker: "ubuntu:latest"
        returnCodes: 10
    }
}


task ReturnCodeString {
    meta {
        volatile: true
    }

    command <<<
        exit 500
    >>>

    runtime {
        docker: "ubuntu:latest"
        returnCodes: "*"
    }
}
