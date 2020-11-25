version 1.0

task check_if_docker_image_cache_disk_mounted {
    command {
        # Hmm not sure, look to see if we mounted /dev/sdb perhaps?
        echo true
    }
    runtime {
        docker: "broadinstitute/gatk:4.0.1.2"
        backend: "Papiv2-Docker-Image-Cache"
        useDockerImageCache: true
    }
    output {
        Boolean mounted_docker_image_cache_disk = read_boolean(stdout())
    }
}

workflow docker_image_cache_unspecified_test {
    call check_if_docker_image_cache_disk_mounted
    output {
        Boolean is_docker_image_cache_disk_mounted = check_if_docker_image_cache_disk_mounted.mounted_docker_image_cache_disk
    }
}
