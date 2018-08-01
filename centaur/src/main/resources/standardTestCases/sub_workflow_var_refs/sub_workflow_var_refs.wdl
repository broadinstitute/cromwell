import "sub_workflow_var_refs_import.wdl" as subworkflow


workflow top_level_workflow {

  Array[String] is = ["fee", "fi", "fo", "fum"]
  Array[String] pieces = ["hello", "lowly", "subject!"]

  scatter (i in is) {
     call subworkflow.subhello as subhello {
      input: greeting_pieces = pieces
    }
  }

  output {
    Array[Int] sal_len_inner = subhello.sal_len
    Int sal_len_outer = length(subhello.hello_out[0])
  }
}
