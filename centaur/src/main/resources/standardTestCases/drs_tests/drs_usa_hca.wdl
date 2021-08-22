version 1.0

workflow drs_usa_hca {
    input {
        File file
    }

    call localize_drs_with_usa {
        input:
            file = file
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
        disks: "local-disk 100 HDD"
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
