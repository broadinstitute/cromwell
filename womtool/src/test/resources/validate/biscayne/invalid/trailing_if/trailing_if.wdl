version development-1.1

workflow trailing_if {
  input {
    Boolean is_exome
  }

  if (is_exome != ) {
    Int is_exome_int = 1
  }
}
