version 1.0

task no_cache {
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

workflow no_cache_delete {
    call no_cache {
        input:
    }

    output {
        File keep = no_cache.keep
    }
}
