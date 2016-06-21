task simpleNeasy {
    command {
       echo "This test is easy, breezy"
       sleep 10
    }
    output {
        String out = read_string(stdout())
    }
    runtime {
        docker: "ubuntu:latest"
    }
}
workflow monitoring {
  call simpleNeasy
  output {
    simpleNeasy.out
  }
}
