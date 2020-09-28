version 1.0

workflow drs_usa_jdr {
    call localize_jdr_drs_with_usa

    output {
        String path1 = localize_jdr_drs_with_usa.path1
        String hash1 = localize_jdr_drs_with_usa.hash1
        Float size1 = localize_jdr_drs_with_usa.size1
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
