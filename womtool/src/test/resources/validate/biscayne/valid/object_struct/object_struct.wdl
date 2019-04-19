version development

struct Object {
  String field
}

workflow object_struct {
  input {
    Object foo
  }

  output {
    Object foo_out = foo
  }
}
