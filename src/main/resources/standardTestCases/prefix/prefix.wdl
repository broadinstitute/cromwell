task pfx {
  Array[String] kvs = ["k1=v1", "k2=v2", "k3=v3"]
  Array[String] prefixed = prefix("-e ", kvs)
  command {
    echo "${sep=' ' prefixed}"
  }
  output {
    String out = read_string(stdout())
  }
}

workflow prefix {
  call pfx
  output {
    String out = pfx.out
  }
}
