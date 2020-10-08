version 1.0

workflow drs_usa_hca {
    input {
        File file1
        File file2
        File file3
        File file4
        File file5
    }

    call localize_drs_with_usa {
        input:
            file1 = file1,
            file2 = file2,
            file3 = file3,
            file4 = file4,
            file5 = file5
    }

    call skip_localize_drs_with_usa {
        input:
            file1 = file1,
            file2 = file2,
            file3 = file3,
            file4 = file4,
            file5 = file5
    }

    output {
        String path1 = localize_drs_with_usa.path1
        String path2 = localize_drs_with_usa.path2
        String path3 = localize_drs_with_usa.path3
        String path4 = localize_drs_with_usa.path4
        String path5 = localize_drs_with_usa.path5
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
        String cloud1 = skip_localize_drs_with_usa.path1
        String cloud2 = skip_localize_drs_with_usa.path2
        String cloud3 = skip_localize_drs_with_usa.path3
        String cloud4 = skip_localize_drs_with_usa.path4
        String cloud5 = skip_localize_drs_with_usa.path5
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
       echo ~{file1} > path1
       echo ~{file2} > path2
       echo ~{file3} > path3
       echo ~{file4} > path4
       echo ~{file5} > path5
       md5sum ~{file1} | cut -c1-32 > hash1
       md5sum ~{file2} | cut -c1-32 > hash2
       md5sum ~{file3} | cut -c1-32 > hash3
       md5sum ~{file4} | cut -c1-32 > hash4
       md5sum ~{file5} | cut -c1-32 > hash5
    >>>

    output {
        String path1 = read_string("path1")
        String path2 = read_string("path2")
        String path3 = read_string("path3")
        String path4 = read_string("path4")
        String path5 = read_string("path5")
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

task skip_localize_drs_with_usa {
    input {
        File file1
        File file2
        File file3
        File file4
        File file5
    }

    parameter_meta {
        file1: { localization_optional: true }
        file2: { localization_optional: true }
        file3: { localization_optional: true }
        file4: { localization_optional: true }
        file5: { localization_optional: true }
    }

    command <<<
        echo ~{file1} > path1
        echo ~{file2} > path2
        echo ~{file3} > path3
        echo ~{file4} > path4
        echo ~{file5} > path5
    >>>

    output {
        String path1 = read_string("path1")
        String path2 = read_string("path2")
        String path3 = read_string("path3")
        String path4 = read_string("path4")
        String path5 = read_string("path5")
    }

    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa"
    }
}
