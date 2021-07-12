task hello {
    command {
        echo "Hello Google Artifact Registry!"
    }
    output {
        String salutation = read_string(stdout())
    }
    runtime {
        docker: "us-central1-docker.pkg.dev/broad-dsde-cromwell-dev/bt-335/ubuntu:bt-335"
    }
}

workflow wf_hello {
    call hello
    output {
        hello.salutation
    }
}
