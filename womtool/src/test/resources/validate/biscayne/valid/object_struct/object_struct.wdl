version development-1.1

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
