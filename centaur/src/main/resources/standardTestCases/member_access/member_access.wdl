task echo_str {
  Array[String] s
  command { 
    echo ${s[0]} ${s[1]} ${s[2]}
  }
  output {
    Pair[String, Int] left = (read_string(stdout()), 27)
  }
  runtime { 
   docker: "ubuntu:latest"
  }
}

workflow test {
  Map[String, Int] m = {"a": 0, "b": 1, "c": 2}
  Array[String] a = ["foo", "bar", "baz"]
  Pair[String, String] p = ("left", "right")
  Pair[Pair[Int, Int], Pair[Int, Int]] triple = ((1, 2), (3, 4))

  call echo_str {
    input: s = [a[1], m["c"], p.left]
  }

  call echo_str as echo_str_2 {
    input: s = [ echo_str.left.left, echo_str.left.right, echo_str.left.left ]
  }

  output {
    # This gets the 'left' output from task echo_str, and then a member lookup for the pair's 'left' element
    String left_out = echo_str.left.left + " " + echo_str_2.left.right

    # Oh hey, nested Pair lookups also works now! Magic!
    Int triple_left = triple.left.left + triple.right.left
  }
}
