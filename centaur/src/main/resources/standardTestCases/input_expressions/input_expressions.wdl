workflow input_expressions {
    Float size256k = size("gs://cloud-cromwell-dev/file_256k", "KiB")
    output {
        Float wf_out = size256k
    }
}
