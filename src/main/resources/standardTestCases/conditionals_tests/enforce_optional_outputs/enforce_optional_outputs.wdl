task double {
  Int x
  command {
    echo $(( ${x} * 2 ))

    sleep 2
  }
  output {
    Int doubled = read_int(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow enforce_optional_outputs {

  if (true) {
    call double as double_i { input: x = 25 }
  }

  output {
    # Error: Should be declared as an optional array:
    Array[Int] oopsNotOptionalArray = [ double_i.doubled, double_i.doubled ]
  }

}
