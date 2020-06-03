version development

struct OptionalIntAndInt {
  Int? optional
  Int non_optional
}

workflow none_literal {
  input {
    Int? int_in = None
  }

  call none_in_command

  output {
    Int? out = int_in
    Int? out2 = None
    Int out3 = select_first([out, out2, None, 3])

    Int? none_from_struct = none_in_command.none_in_struct.optional
    String out_none = none_in_command.s
  }

}

task none_in_command {
  command <<<
    echo ~{if (defined(None)) then "Some" else "None" }
    echo 4 > 4.txt
  >>>
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    String s = read_string(stdout())
    OptionalIntAndInt none_in_struct = object { optional: None, non_optional: read_int("4.txt") }
  }
}
