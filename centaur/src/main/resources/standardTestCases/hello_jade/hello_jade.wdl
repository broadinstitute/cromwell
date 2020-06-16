version 1.0

workflow hello_jade {
    # Uses test file https://jade.datarepo-dev.broadinstitute.org/snapshots/details/93dc1e76-8f1c-4949-8f9b-07a087f3ce7b
    # The file takes 10 minutes to localize. ASAP we should upload a smaller file and update the expected .test results.
    # `md5sum` of the 15gb file also takes too long. Uncomment the `md5sum` below when we switch files.
    File jade_file = "drs://jade.datarepo-dev.broadinstitute.org/v1_93dc1e76-8f1c-4949-8f9b-07a087f3ce7b_8b07563a-542f-4b5c-9e00-e8fe6b1861de"
    call go_jade { input: jade_file = jade_file }
    output { Array[String] out = go_jade.out }
}

task go_jade {
    input { File jade_file }
    output { Array[String] out = read_lines("out.txt") }
    command <<<
        # md5sum '~{jade_file}' | awk '{print $1}' >> out.txt
        du '~{jade_file}' | awk '{print $1}' >> out.txt
        ls '~{jade_file}' >> out.txt
    >>>
    runtime {
        docker: "ubuntu"
        backend: "papi-v2-usa" # As of 2020-06 this is one of the only centaur backends with dos/drs
        disks: "local-disk " + ceil(size(jade_file, "GB") + 10) + " SSD"
    }
}
