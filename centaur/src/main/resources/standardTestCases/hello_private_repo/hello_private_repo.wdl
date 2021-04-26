version 1.0

task hello {
    input {
        String addressee
    }
    command {
        echo "Hello ~{addressee}!"
    }
    output {
        String salutation = read_string(stdout())
    }
    runtime {
        backend: "LocalDockerSecure"
        docker: "broadinstitute/cromwell-docker-test:private-repo"
    }
}

workflow hello_private_repo {
    input {
        String addressee
    }
    call hello { input: addressee = addressee }
    output {
        String salutation = hello.salutation
    }
}
