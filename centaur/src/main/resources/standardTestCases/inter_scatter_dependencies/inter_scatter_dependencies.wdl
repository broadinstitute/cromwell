task mirror {
  Array[String] a
  command {
    echo ${sep="_" a}
  }
  output {
    String out = read_string(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task mirror_array {
  Array[String] a
  command {
    for i in ${sep=" "a}
    do
      echo $i
    done
  }
  output {
    Array[String] out = read_lines(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow inter_scatter_dependencies {
  String hello = "hello"

  Array[Int] int_array = range(5)

  scatter(a in int_array) {
    call mirror as mirror1 { input: a = [hello, "world"] }
    call mirror as mirror3 { input: a = mirror2.out }
  }

  call mirror_array as mirror2 { input: a = mirror1.out }

  output {
    Array[String] result = mirror3.out
  }
}
