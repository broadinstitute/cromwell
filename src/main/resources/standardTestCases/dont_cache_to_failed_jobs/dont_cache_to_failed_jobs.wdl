task failing_task {
    command {
        # C'est nul !
        exit 1
    }
    runtime {
        docker: "ubuntu@sha256:a0ee7647e24c8494f1cf6b94f1a3cd127f423268293c25d924fbe18fd82db5a4"
    }
    output {
        Boolean done = true
    }
}

workflow dont_cache_to_failed_jobs {
  call failing_task
}
