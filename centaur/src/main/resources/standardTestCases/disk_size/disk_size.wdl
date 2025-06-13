version 1.0

workflow disk_size {
    call ds_optional as optional_provided {
        input:
            disk_size = 20
    }

    call ds_optional as optional_omitted

    call ds_required as required_provided {
        input:
            disk_size = 40
    }

    call ds_required as required_omitted
}

task ds_optional {
    input {
        Int? disk_size = 60
    }

    command {
        echo "Disk size is ${disk_size} GB"
    }

    output {
        String message = read_string(stdout())
    }

    runtime {
        docker: "ubuntu:latest"
        bootDiskSizeGb: disk_size
    }
}

task ds_required {
    input {
        Int disk_size = 80
    }

    command {
        echo "Disk size is ${disk_size} GB"
    }

    output {
        String message = read_string(stdout())
    }

    runtime {
        docker: "ubuntu:latest"
        bootDiskSizeGb: disk_size
    }
}
