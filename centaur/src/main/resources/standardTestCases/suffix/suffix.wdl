task sfx {
  Array[String] filenames = ["logs", "cmd", "output"]
  Array[String] suffixed = suffix(".txt", filenames)
  command {
    echo "${sep=' ' suffixed}"
  }
  output {
    String out = read_string(stdout())
  }
  runtime { docker: "ubuntu:latest" }
}

workflow suffix {
  call sfx
  output {
    String out = sfx.out
  }
}
