version 1.0

struct PutMeInAPair {
  Int value
}

workflow member_access {
  Map[String, Int] m = {"a": 0, "b": 1, "c": 2}
  Array[String] a = ["foo", "bar", "baz"]
  Pair[String, String] p = ("left", "right")
  Pair[Pair[Int, Int], Pair[Int, Int]] triple = ((1, 2), (3, 4))
  Pair[Pair[Array[Int], Int], Int] array_in_pair = (([1,2,3], 4), 5)
  PutMeInAPair pmiap1 = object { value: 1111 }
  PutMeInAPair pmiap2 = object { value: 222 }
  PutMeInAPair pmiap3 = object { value: 333 }
  Pair[Pair[Array[PutMeInAPair], Int], Int] struct_in_pair = (([pmiap1, pmiap2, pmiap3], 4), 5)

  call echo_str {
    input: s = [a[1], m["c"], p.left]
  }

  call echo_str as echo_str_2 {
    input: s = [ echo_str.left.left, echo_str.left.right, echo_str.left.left ]
  }

  output {
    String left_out = echo_str.left.left + " " + echo_str_2.left.right
    Int triple_left = triple.left.left + triple.right.left
    Int array_in_pair_lookup = array_in_pair.left.left[1]
    Int struct_in_pair_lookup = struct_in_pair.left.left[1].value
  }
}

task echo_str {
  input {
    Array[String] s
  }
  command { 
    echo ~{s[0]} ~{s[1]} ~{s[2]}
  }
  output {
    Pair[String, Int] left = (read_string(stdout()), 27)
  }
  runtime { 
   docker: "ubuntu:latest"
  }
}
