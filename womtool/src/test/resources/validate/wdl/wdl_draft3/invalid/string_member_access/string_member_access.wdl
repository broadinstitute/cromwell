version 1.0

workflow string_member_access {
  String s = "hello"

  output {
    Int x = "~{s}".tbi
  }
}
