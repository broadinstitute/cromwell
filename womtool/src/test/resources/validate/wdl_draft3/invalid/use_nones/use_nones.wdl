version 1.0

workflow use_nones {
  input {
    Int? int_in = None
  }

  output {
    Int? out = int_in
    Int? out2 = None
  }

}
