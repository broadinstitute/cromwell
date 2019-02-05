version 1.0

task mkdir {
    input {
        String dir_name
        Boolean empty
    }
    String touch1 = if empty then "" else "touch '~{dir_name}'/file1"
    String touch2 = if empty then "" else "touch '~{dir_name}'/file2"
    command {
        mkdir '~{dir_name}'
        ~{touch1}
        ~{touch2}
    }
    output {
        File empty_dir = dir_name
    }
    runtime {
        docker: "ubuntu"
    }
}
workflow empty_dir_workflow {
    call mkdir as mkdir_normal_empty { input: dir_name = 'empty_normal_dir', empty = true }
    call mkdir as mkdir_dollar_empty { input: dir_name = 'empty_$dollar_dir', empty = true }
    call mkdir as mkdir_space_empty { input: dir_name = 'empty_space dir', empty = true }
    call mkdir as mkdir_normal_nonempty { input: dir_name = 'nonempty_normal_dir', empty = false }
    call mkdir as mkdir_dollar_nonempty { input: dir_name = 'nonempty_$dollar_dir', empty = false }
    call mkdir as mkdir_space_nonempty { input: dir_name = 'nonempty_space dir', empty = false }
}