version 1.0

task generate {
    input {
    }

    command <<<
        echo "I wanna be already deleted" > delete.txt
    >>>

    runtime {
        docker: "ubuntu"
    }

    output {
        File delete = "delete.txt"
    }
}

task delete {
    input {
        String path
    }

    command <<<
        gsutil rm '~{path}'
    >>>

    runtime {
        docker: "google/cloud-sdk:latest"
    }

    output {
    }
}

workflow missing_delete {
    call generate {
        input:
    }
    call delete {
        input:
            path = generate.delete
    }

    output {
    }
}
