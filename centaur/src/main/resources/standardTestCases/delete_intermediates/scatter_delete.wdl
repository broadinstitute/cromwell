version 1.0

task scatterer {
    input {
        Int index
    }

    command <<<
        echo "I wanna be deleted ~{index}" > delete.txt
        echo "I think I'll go for a walk ~{index}" > keep.txt
    >>>

    runtime {
        docker: "ubuntu"
    }

    output {
        File delete = "delete.txt"
        File keep = "keep.txt"
    }
}

workflow scatter_delete {
    scatter(index in range(3)) {
        call scatterer {
            input: index = index
        }
    }

    output {
        Array[File] kept = scatterer.keep
    }
}
