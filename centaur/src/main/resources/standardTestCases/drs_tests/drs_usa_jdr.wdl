version 1.0

workflow drs_usa_jdr {
    input {
        File file1
    }

    call localize_jdr_drs_with_usa {
        input:
            file1 = file1
    }

    call skip_localize_jdr_drs_with_usa {
        input:
            file1 = file1
    }

    call read_drs_with_usa {
        input:
            file1 = file1
    }

    output {
        String path1 = localize_jdr_drs_with_usa.path1
        String hash1 = localize_jdr_drs_with_usa.hash1
        Float size1 = localize_jdr_drs_with_usa.size1
        String cloud1 = skip_localize_jdr_drs_with_usa.path1
        Map[String, String] map1 = read_drs_with_usa.map1
    }
}

task localize_jdr_drs_with_usa {
    input {
       File file1
    }

    command <<<
       echo ~{file1} > path1
       md5sum ~{file1} | cut -c1-32 > hash1
    >>>

    output {
        String path1 = read_string("path1")
        String hash1 = read_string("hash1")
        Float size1 = size(file1)
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }
}

task skip_localize_jdr_drs_with_usa {
    input {
       File file1
    }

    parameter_meta {
        file1: { localization_optional: true }
    }

    command <<<
       echo ~{file1} > path1
    >>>

    output {
        String path1 = read_string("path1")
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }
}

task read_drs_with_usa {
    input {
        File file1
    }

    command <<<
        echo file is read by the engine
    >>>

    output {
        Map[String, String] map1 = read_json(file1)
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }
}
