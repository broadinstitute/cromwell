task goodbye {
  String addressee
  command {
    echo "Goodbye ${addressee}!"
  }
  output {
    String farewell = read_string(stdout())
  }
}

workflow local_backend {
  call goodbye
  output {
     goodbye.farewell
  }
}
