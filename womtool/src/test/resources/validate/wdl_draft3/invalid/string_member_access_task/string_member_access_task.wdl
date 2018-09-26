version 1.0

task string_member_access_output {
  String s = "hello"
  command {}
  output {
    String out = "~{s}".tbi
  }
}

task string_member_access_output_runtime {
  String s = "hello"
  command {}
  runtime {
    docker: "~{s}".latest
  }
  output {
    String out = "~{s}"
  }
}

task x
{
  input {
    String x = "abc"
  }
  command {}
  output {
    File y = "${x}".xyz
  }
}
