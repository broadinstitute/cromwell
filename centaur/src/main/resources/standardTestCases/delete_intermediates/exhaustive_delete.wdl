version 1.0

task exhaustive {
    input {
    }

    command <<<
        echo "All other files will be listed in outputExpectations, like tears in rain" > delete.txt
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

workflow exhaustive_delete {
    call exhaustive {
        input:
    }

    output {
        File keep = exhaustive.keep
    }
}
