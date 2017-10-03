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
    Int single = 25
    call double as double_i { input: x = single }
  }
  if (false) {
    call double as double_j { input: x = 26 }
  }

  output {
    Array[Int] oopsNotOptionalArray = [ double_i.doubled, double_j.doubled ]
  }

}
