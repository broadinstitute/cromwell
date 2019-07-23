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
    # For all below 5 HCA uuids, Martha does not return a service account
    File file1 = "dos://service.staging.explore.data.humancellatlas.org/033c9840-c5cd-438b-b0e4-8e4cd8fc8dc6?version=2019-07-04T104122.106166Z"
    File file2 = "dos://service.staging.explore.data.humancellatlas.org/4defa7b0-46c2-4053-8e99-b827eed1bc96?version=2019-07-04T104122.100969Z"
    File file3 = "dos://service.staging.explore.data.humancellatlas.org/de5dcfc1-5aea-41ba-a7ae-e72c416cb450?version=2019-07-04T104122.092788Z"
    File file4 = "dos://service.staging.explore.data.humancellatlas.org/16dea2c5-e2bd-45bc-b2fd-fcac0daafc48?version=2019-07-04T104122.060634Z"
    File file5 = "dos://service.dev.explore.data.humancellatlas.org/7c800467-9143-402f-b965-4e7cad75c1e6?version=2019-05-26T130511.722646Z"

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
