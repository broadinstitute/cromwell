task bigtask {
    File in_file

    command {
        du -h ${in_file}
    }

    runtime {
        docker: "ubuntu:latest"
        disks: "local-disk 100 SSD"
    }

    output {
        File out = in_file
    }
}

task otherbigtask {
    File in_file

    command {
        du -h ${in_file}
    }

    runtime {
        docker: "ubuntu:latest"
        disks: "local-disk 100 SSD"
    }
}

workflow large_input {
    call bigtask
    call otherbigtask {input: in_file = bigtask.out}
}
