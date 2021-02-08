version 1.0

workflow drs_usa_jdr {
    input {
        File file1
        File file2
    }

    call localize_jdr_drs_with_usa {
        input:
            file1 = file1,
            file2 = file2
    }

    call skip_localize_jdr_drs_with_usa {
        input:
            file1 = file1,
            file2 = file2
    }

    call read_drs_with_usa {
        input:
            file1 = file1,
            file2 = file2
    }

    output {
        String path1 = localize_jdr_drs_with_usa.path1
        String path2 = localize_jdr_drs_with_usa.path2
        String hash1 = localize_jdr_drs_with_usa.hash1
        String hash2 = localize_jdr_drs_with_usa.hash2
        Float size1 = localize_jdr_drs_with_usa.size1
        Float size2 = localize_jdr_drs_with_usa.size2
        String cloud1 = skip_localize_jdr_drs_with_usa.path1
        String cloud2 = skip_localize_jdr_drs_with_usa.path2
        Map[String, String] map1 = read_drs_with_usa.map1
        Map[String, String] map2 = read_drs_with_usa.map2
    }
}

task localize_jdr_drs_with_usa {
    input {
       File file1
       File file2
    }

    command <<<
       echo ~{file1} > path1
       echo ~{file2} > path2
       md5sum ~{file1} | cut -c1-32 > hash1
       md5sum ~{file2} | cut -c1-32 > hash2
    >>>

    output {
        String path1 = read_string("path1")
        String path2 = read_string("path2")
        String hash1 = read_string("hash1")
        String hash2 = read_string("hash2")
        Float size1 = size(file1)
        Float size2 = size(file2)
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }
}

task skip_localize_jdr_drs_with_usa {
    input {
       File file1
       File file2
    }

    parameter_meta {
        file1: { localization_optional: true }
        file2: { localization_optional: true }
    }

    command <<<
       echo ~{file1} > path1
       echo ~{file2} > path2
    >>>

    output {
        String path1 = read_string("path1")
        String path2 = read_string("path2")
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }
}

task read_drs_with_usa {
    input {
        File file1
        File file2
    }

    command <<<
        echo file is read by the engine
    >>>

    output {
        Map[String, String] map1 = read_json(file1)
        Map[String, String] map2 = read_json(file2)
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }
}
