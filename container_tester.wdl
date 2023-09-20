version 1.0

workflow container_tester {

  input {
    String image_name
    String uuid
  }

  call container_tester_task {
    input:
      image_name = image_name,
      uuid = uuid
  }

  output {
    String result = container_tester_task.out
  }
}

task container_tester_task {
  input {
    String uuid
    String image_name
  }

  command {
    echo "I am doing absolutely nothing with image ~{image_name} and UUID ~{uuid}"
  }

  output {
    String out = read_string(stdout())
  }

  runtime {
    docker: docker_image
  }
}
