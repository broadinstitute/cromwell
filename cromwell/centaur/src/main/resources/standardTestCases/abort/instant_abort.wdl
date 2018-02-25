task aborted {
    command {
        echo "Instant abort"
    }
    
    runtime {
        docker: "ubuntu:latest"
    }
}

workflow instant_abort {
    call aborted
}
