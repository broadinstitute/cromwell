task echo_str {
  String s
  command { echo ${s} }
  output { String o = read_string(stdout()) }
  runtime { 
   docker: "ubuntu:latest"
  }
}

task echo_int {
  Int i
  command { echo ${i} }
  output { Int o = read_int(stdout()) }
  runtime {
   docker: "ubuntu:latest"
  }
}

workflow test {
  Map[String, Int] m = {"a": 0, "b": 1, "c": 2}
  Array[String] a = ["foo", "bar", "baz"]
  call echo_str {
    input: s = a[1]
  }
  call echo_int {
    input: i = m["c"]*100
  }
}
