version 1.0

workflow drs_usa_hca {
    call localize_drs_with_usa

    output {
        String hash1 = localize_drs_with_usa.hash1
        String hash2 = localize_drs_with_usa.hash2
        String hash3 = localize_drs_with_usa.hash3
        String hash4 = localize_drs_with_usa.hash4
        String hash5 = localize_drs_with_usa.hash5
        Float size1 = localize_drs_with_usa.size1
        Float size2 = localize_drs_with_usa.size2
        Float size3 = localize_drs_with_usa.size3
        Float size4 = localize_drs_with_usa.size4
        Float size5 = localize_drs_with_usa.size5
    }
}

task localize_drs_with_usa {
    input {
       File file1
       File file2
       File file3
       File file4
       File file5
    }

    command <<<
       md5sum ~{file1} | cut -c1-32 > hash1
       md5sum ~{file2} | cut -c1-32 > hash2
       md5sum ~{file3} | cut -c1-32 > hash3
       md5sum ~{file4} | cut -c1-32 > hash4
       md5sum ~{file5} | cut -c1-32 > hash5
    >>>

    output {
        String hash1 = read_string("hash1")
        String hash2 = read_string("hash2")
        String hash3 = read_string("hash3")
        String hash4 = read_string("hash4")
        String hash5 = read_string("hash5")
        Float size1 = size(file1)
        Float size2 = size(file2)
        Float size3 = size(file3)
        Float size4 = size(file4)
        Float size5 = size(file5)
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }
}
