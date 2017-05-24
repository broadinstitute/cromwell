task expression_locations {
  Int i
  Int j
  # In an input declaration
  Int max = if i > j then i else j

  # In a command:
  command {
    echo ${max} > max.txt
    echo ${if i < j then i else j} > min.txt
  }
  output {
    # In an output declaration:
    String maxOrMin = if (i % 2 == 0) then read_string("max.txt") else read_string("min.txt")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task only_evaluate_correct_side {
  Boolean b = true
  command {
    echo t > t
    echo f > f
  }
  output {
    # This'll fail if we try to evaluate the wrong side!
    String only_lhs_works = if b then read_string("t") else read_int("t")
    String only_rhs_works = if !b then read_int("f") else read_string("f")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow if_then_else_expressions {
  Int five = if "A" == "A" then 5 else 6
  call expression_locations as fiveVsTen { input: i = five, j = 10 }
  call expression_locations as sixVsTen { input: i = if "A" == "B" then 5 else 6, j = 10 }
  call only_evaluate_correct_side

  output {
    String a = fiveVsTen.maxOrMin
    String b = sixVsTen.maxOrMin
    String c = if a == b then only_evaluate_correct_side.only_lhs_works else only_evaluate_correct_side.only_rhs_works
    String d = if a != b then only_evaluate_correct_side.only_lhs_works else only_evaluate_correct_side.only_rhs_works
  }
}
