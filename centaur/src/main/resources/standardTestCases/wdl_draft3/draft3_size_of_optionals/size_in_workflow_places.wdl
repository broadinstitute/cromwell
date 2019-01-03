version 1.0

workflow size_in_workflow_places {
  input {
    File? input_file_1
    File? input_file_2
    Array[File?] input_files

    Int file_size_1_in_inputs = round(size(input_file_1))
    Int file_size_2_in_inputs = round(size(input_file_2))
    Int array_size_in_inputs = round(size(input_files))
  }

  Int file_size_1_in_declaration = round(size(input_file_1))
  Int file_size_2_in_declaration = round(size(input_file_2))
  Int array_size_in_declaration = round(size(input_files))

  output {
    Int file_size_1_in_inputs_out = file_size_1_in_inputs
    Int file_size_2_in_inputs_out = file_size_2_in_inputs
    Int array_size_in_inputs_out = array_size_in_inputs

    Int file_size_1_in_declaration_out = file_size_1_in_declaration
    Int file_size_2_in_declaration_out = file_size_2_in_declaration
    Int array_size_in_declaration_out = array_size_in_declaration

    Int file_size_1_in_outputs = round(size(input_file_1))
    Int file_size_2_in_outputs = round(size(input_file_2))
    Int array_size_in_outputs = round(size(input_files))
  }
}
