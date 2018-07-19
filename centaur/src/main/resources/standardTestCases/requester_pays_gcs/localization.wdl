version 1.0

workflow requester_pays_localization {
    File input_file = "gs://cromwell_bucket_with_requester_pays/lorem ipsum.txt"
    call localize { input: f = input_file }
    output {
        String content = localize.o
        Array[String] globbed = localize.g
    }
}

task localize {
    input {
        File f
    }
    command {
        cat "~{f}" > out.txt
    }
    runtime {
        docker: "ubuntu"
    }
    output {
        String o = read_string(stdout())
        Array[String] g = glob("*.txt")
    }
}
