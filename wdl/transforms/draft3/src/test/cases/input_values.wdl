version 1.0

workflow input_values {
  input {
    # All the primitive types:
    Int i = 5
    String s = "hello"
    String placeholder = "h${s}o"
    String placeholder2 = 'h~{s}o'
    Float f = 5.5
    Boolean b = true
  }
}
