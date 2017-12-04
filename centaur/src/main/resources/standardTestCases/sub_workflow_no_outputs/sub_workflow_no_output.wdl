import "sub_workflow_no_output_block_import.wdl" as no_output_block
import "sub_workflow_no_outputs_in_block_import.wdl" as no_outputs_in_block
import "sub_workflow_one_output_import.wdl" as one_output

workflow top_level_workflow {

  Array[Int] is = range(4)
  Array[String] pieces = ["hello", "lowly", "subject!"]

  call no_output_block.subhello as nob { input: greeting_pieces = pieces }
  call no_outputs_in_block.subhello as noib { input: greeting_pieces = pieces }
  call one_output.subhello as oo { input: greeting_pieces = pieces }

  # Despite no outputs here, we DO expect an output to appear from one_output
  output { }
}
