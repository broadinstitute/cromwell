import "sub_workflow_var_refs_import.wdl" as subworkflow


workflow top_level_workflow {

  Array[String] is = ["fee", "fi", "fo", "fum"]
  Array[String] pieces = ["hello", "lowly", "subject!"]

  scatter (i in is) {
    call subworkflow.subhello as subhello_call_alias {
      input: greeting_pieces = pieces
    }
  }

  output {
    Array[Int] sal_len_inner = subhello_call_alias.sal_len
    Int sal_len_outer = length(subhello_call_alias.hello_out[0])
  }
}
