version draft-3

workflow my_workflow {
  input{
    Int i
  }
  output {
    Int o = i + 100
  }
}
