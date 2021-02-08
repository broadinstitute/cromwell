version 1.0

task hello {
    input {
    }

    command <<<
        echo "I wanna be deleted" > delete.txt
        echo "I think I'll go for a walk" > keep.txt
    >>>

    runtime {
        docker: "ubuntu"
    }

    output {
        File delete = "delete.txt"
        File keep = "keep.txt"
    }
}

workflow hello_delete {
    call hello {
        input:
    }

    output {
        File keep = hello.keep
    }
}
