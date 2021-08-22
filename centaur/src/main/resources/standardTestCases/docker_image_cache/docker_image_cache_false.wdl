version 1.0

task check_if_docker_image_cache_disk_mounted {
  String google_docker_cache_disk_name = "google-docker-cache"
  command {
    instance_metadata_disks=$(curl "http://metadata.google.internal/computeMetadata/v1/instance/disks/?recursive=true" -H "Metadata-Flavor: Google")
    if [[ "$instance_metadata_disks" == *~{google_docker_cache_disk_name}* ]]; then
      echo "true"
    else
      echo "false"
    fi
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.3.0"
    backend: "Papiv2-Docker-Image-Cache"
    useDockerImageCache: false
  }
  meta {
    volatile: true
  }
  output {
    Boolean mounted_docker_image_cache_disk = read_boolean(stdout())
  }
}

workflow docker_image_cache_false_test {
  call check_if_docker_image_cache_disk_mounted
  output {
    Boolean is_docker_image_cache_disk_mounted = check_if_docker_image_cache_disk_mounted.mounted_docker_image_cache_disk
  }
}
