task bad {
  command {
    echo "hello" > a
    sleep 2
  }

  runtime {
    docker: "ubuntu:latest"
  }

  output {
    # Oops! we made a spelling mistake in our WDL!
    File a = "b"
  }
}

workflow badExample {
  call bad
}
