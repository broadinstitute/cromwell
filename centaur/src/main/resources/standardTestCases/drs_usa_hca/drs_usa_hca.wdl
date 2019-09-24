version 1.0

workflow drs_usa_hca {
    call localize_drs_with_usa

    output {
        String hash1 = localize_drs_with_usa.hash1
        String hash2 = localize_drs_with_usa.hash2
        String hash3 = localize_drs_with_usa.hash3
        Float size1 = localize_drs_with_usa.size1
        Float size2 = localize_drs_with_usa.size2
        Float size3 = localize_drs_with_usa.size3
    }
}

task localize_drs_with_usa {
    input {
       File file1
       File file2
       File file3
    }

    command <<<
       md5sum ~{file1} | cut -c1-32 > hash1
       md5sum ~{file2} | cut -c1-32 > hash2
       md5sum ~{file3} | cut -c1-32 > hash3
    >>>

    output {
        String hash1 = read_string("hash1")
        String hash2 = read_string("hash2")
        String hash3 = read_string("hash3")
        Float size1 = size(file1)
        Float size2 = size(file2)
        Float size3 = size(file3)
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }
}
