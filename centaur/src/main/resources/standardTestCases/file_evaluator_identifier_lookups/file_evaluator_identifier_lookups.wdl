version 1.0

workflow file_evaluator_identifier_lookups {
  call identifier
  call identifier_member_access
  call index_access
}

task identifier {
  String testFile = "test.txt"

  command {
    echo OK > ~{testFile}
  }

  output {
    String id = read_string(testFile)
  }

  runtime {
    docker: "debian:stable-slim"
  }
}

task identifier_member_access {
  Pair[String, String] pair = ("foo.txt", "bar.txt")

  command {
    echo LEFT > ~{pair.left}
    echo RIGHT > ~{pair.right}
  }

  output {
    String left = read_string(pair.left)
    String right = read_string(pair.right)
  }

  runtime {
    docker: "debian:stable-slim"
  }
}

task index_access {
  Array[String] array = ["foo.txt", "bar.txt"]

  command {
    echo ZERO > ~{array[0]}
    echo ONE > ~{array[1]}
  }

  output {
    String zero = read_string(array[0])
    String one = read_string(array[1])
  }

  runtime {
    docker: "debian:stable-slim"
  }
}
