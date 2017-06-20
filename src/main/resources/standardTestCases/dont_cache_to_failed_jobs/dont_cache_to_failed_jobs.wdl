task failing_task {
    command {
        # C'est nul !
        exit 1
    }
    runtime {
        docker: "ubuntu:latest"
    }
    output {
        Boolean done = true
    }
}

workflow dont_cache_to_failed_jobs {
  call failing_task
}
