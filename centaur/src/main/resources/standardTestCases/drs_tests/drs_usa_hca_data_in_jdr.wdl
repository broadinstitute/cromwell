version 1.0

workflow drs_usa_hca_data_in_jdr {
    input {
        File file
        Int disk_size_required_for_file_gb
    }

    call localize_drs_with_usa {
        input:
            file = file,
            disk_size_required_for_file_gb = disk_size_required_for_file_gb
    }

    call skip_localize_drs_with_usa { input: file = file }

    output {
        String path = localize_drs_with_usa.path
        String hash = localize_drs_with_usa.hash
        Float size = localize_drs_with_usa.size
        String cloud = skip_localize_drs_with_usa.path
    }
}

task localize_drs_with_usa {
    input {
        File file
        Int disk_size_required_for_file_gb
    }

    command <<<
        echo ~{file} > path
        md5sum ~{file} | cut -c1-32 > hash
    >>>

    output {
        String path = read_string("path")
        String hash = read_string("hash")
        Float size = size(file)
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
        disks: "local-disk ~{disk_size_required_for_file_gb} SSD"
    }
}

task skip_localize_drs_with_usa {
    input {
        File file
    }

    parameter_meta {
        file: { localization_optional: true }
    }

    command <<<
        echo ~{file} > path
    >>>

    output {
        String path = read_string("path")
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }
}
