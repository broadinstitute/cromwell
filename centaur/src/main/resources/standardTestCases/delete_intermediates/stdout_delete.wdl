version 1.0

task do_stdout {
    input {
    }

    command <<<
        echo "I wanna be deleted and I'm taking stdout with me" > delete.txt
        echo "I think I'll go for a walk" > keep.txt
    >>>

    runtime {
        docker: "ubuntu"
    }

    output {
        File delete = "delete.txt"
        File stdout = stdout()
        File keep = "keep.txt"
    }
}

workflow stdout_delete {
    call do_stdout {
        input:
    }

    output {
        File keep = do_stdout.keep
    }
}
