version 1.2

task env_variables {
  input {
    env String greeting
    env Array[String] names
    String plain = "not exported"
  }

  env Int count = length(names)

  command <<<
    # Values are exported as environment variables, so they are never interpreted by the shell
    echo "$greeting" > greeting.txt
    echo "$names" > names.txt
    echo "$count" > count.txt
    echo "${plain:-unset}" > plain.txt
  >>>

  output {
    String greeting_out = read_string("greeting.txt")
    String names_out = read_string("names.txt")
    Int count_out = read_int("count.txt")
    String plain_out = read_string("plain.txt")
  }

  requirements {
    container: "ubuntu:latest"
  }
}

workflow env_variables_wf {
  call env_variables {
    greeting = "it's $HOME; echo injected",
    names = ["a", "b"]
  }

  output {
    String greeting = env_variables.greeting_out
    String names = env_variables.names_out
    Int count = env_variables.count_out
    String plain = env_variables.plain_out
  }
}
