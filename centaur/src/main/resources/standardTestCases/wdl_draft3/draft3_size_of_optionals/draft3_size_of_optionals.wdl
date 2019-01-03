version 1.0

import "size_in_workflow_places.wdl" as siwp

workflow draft3_size_of_optionals {
  call make_a_file

  if(false) {
    call make_a_file as dont_make_a_file
  }

  scatter (i in range(2)) {
    if (i % 2 == 0) {
      call make_a_file as scattered_conditional_file_making
    }
    Int made_file_size = round(size(scattered_conditional_file_making.a_file))
  }

  call size_in_task_places { input:
    input_file_1 = make_a_file.a_file,
    input_file_2 = dont_make_a_file.a_file,
    input_files = scattered_conditional_file_making.a_file
  }

  call siwp.size_in_workflow_places { input:
    input_file_1 = make_a_file.a_file,
    input_file_2 = dont_make_a_file.a_file,
    input_files = scattered_conditional_file_making.a_file
  }

  output {
    Int immediate_made_file_size = make_a_file.immediate_file_size
    Int? immediate_not_made_file_size = dont_make_a_file.immediate_file_size

    Array[Int] made_file_sizes = made_file_size
    Int total_made_file_size = round(size(scattered_conditional_file_making.a_file))

    Array[String] lines = size_in_task_places.lines
    Int file_size_1_in_outputs = size_in_task_places.file_size_1_in_outputs
    Int file_size_2_in_outputs = size_in_task_places.file_size_2_in_outputs
    Int array_size_in_outputs = size_in_task_places.array_size_in_outputs

    Int siwp_file_size_1_in_inputs_out = size_in_workflow_places.file_size_1_in_inputs_out
    Int siwp_file_size_2_in_inputs_out = size_in_workflow_places.file_size_2_in_inputs_out
    Int siwp_array_size_in_inputs_out = size_in_workflow_places.array_size_in_inputs_out

    Int siwp_file_size_1_in_declaration_out = size_in_workflow_places.file_size_1_in_declaration_out
    Int siwp_file_size_2_in_declaration_out = size_in_workflow_places.file_size_2_in_declaration_out
    Int siwp_array_size_in_declaration_out = size_in_workflow_places.array_size_in_declaration_out

    Int siwp_file_size_1_in_outputs = size_in_workflow_places.file_size_1_in_outputs
    Int siwp_file_size_2_in_outputs = size_in_workflow_places.file_size_2_in_outputs
    Int siwp_array_size_in_outputs = size_in_workflow_places.array_size_in_outputs
  }
}

task make_a_file {
  command {
    echo "this file is 22 bytes" > a_file
  }
  output {
    File a_file = "a_file"
    Int immediate_file_size = round(size(a_file))
  }
  runtime { docker: "ubuntu:latest" }
}

task size_in_task_places {
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

  command {
    echo ~{file_size_1_in_inputs}
    echo ~{file_size_2_in_inputs}
    echo ~{array_size_in_inputs}
    echo ~{file_size_1_in_declaration}
    echo ~{file_size_2_in_declaration}
    echo ~{array_size_in_declaration}
  }

#  echo ~{round(size(input_file_1))}
#  echo ~{round(size(input_file_2))}
#  echo ~{round(size(input_files))}

  runtime { docker: "ubuntu:latest" }

  output {
    Array[String] lines = read_lines(stdout())
    Int file_size_1_in_outputs = round(size(input_file_1))
    Int file_size_2_in_outputs = round(size(input_file_2))
    Int array_size_in_outputs = round(size(input_files))
  }
}
