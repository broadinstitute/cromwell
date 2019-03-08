version 1.0

task localize_file {
    input {
        File input_file
    }
    command {
        cat "localizing file over 1 GB"
    }
    runtime {
        docker: "ubuntu:latest"
        disks: "local-disk 1 HDD"
    }
    output {
        String out = read_string(stdout())
    }
}

workflow localize_file_larger_than_disk_space {
    File wf_input = "gs://cromwell_test_bucket/file_over_1_gb.txt"

    call localize_file { input: input_file = wf_input }

    output {
        String content = localize_file.out
    }
}
