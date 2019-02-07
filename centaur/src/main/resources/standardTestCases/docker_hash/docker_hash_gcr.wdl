task gcr {
    String image
    command {
        echo "hello"
    }
    runtime {
        docker: image
    }
}

workflow docker_hash_gcr {
    call gcr { input: image = "gcr.io/google-containers/ubuntu-slim:0.14" }
    call gcr as gcr_eu { input: image = "eu.gcr.io/google-containers/ubuntu-slim:0.14" }
    call gcr as gcr_us { input: image = "us.gcr.io/google-containers/ubuntu-slim:0.14" }
    call gcr as gcr_asia { input: image = "asia.gcr.io/google-containers/ubuntu-slim:0.14" }
}
