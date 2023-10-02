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

    # Rather than resolving DRS urls to signed URLs, download files from GCP directly.
    # This invokes a different code path in the DRSLocalizer that uses gsutil rather than getm.
    call read_from_gcp_directly {
        input:
            file1 = "gs://broad-jade-dev-data-bucket/ca8edd48-e954-4c20-b911-b017fedffb67/585f3f19-985f-43b0-ab6a-79fa4c8310fc",
            file2 = "gs://broad-jade-dev-data-bucket/e1941fb9-6537-4e1a-b70d-34352a3a7817/ad783b60-aeba-4055-8f7b-194880f37259/hello_jade_2.json"
    }

    # Verify that the localizer can handle a hybrid of drs and GCS urls.
    call read_from_drs_and_gcp {
        input:
            file1 = "gs://broad-jade-dev-data-bucket/ca8edd48-e954-4c20-b911-b017fedffb67/585f3f19-985f-43b0-ab6a-79fa4c8310fc",
            file2 = file2
    }

    call read_single_file {
        input:
            file1 = file1
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
        Map[String, String] map3 = read_from_gcp_directly.map1
        Map[String, String] map4 = read_from_gcp_directly.map2
        Map[String, String] map5 = read_from_drs_and_gcp.map1
        Map[String, String] map6 = read_from_drs_and_gcp.map2
        Map[String, String] map7 = read_single_file.map1
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
        docker: "ubuntu:latest"
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
        docker: "ubuntu:latest"
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
        docker: "ubuntu:latest"
        backend: "papi-v2-usa"
    }
}

task read_from_gcp_directly {
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
        docker: "ubuntu:latest"
        backend: "papi-v2-usa"
    }
}

task read_from_drs_and_gcp {
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
        docker: "ubuntu:latest"
        backend: "papi-v2-usa"
    }
}

task read_single_file {
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
        docker: "ubuntu:latest"
        backend: "papi-v2-usa"
    }
}
    