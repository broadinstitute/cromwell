version 1.0

task file_to_lines {
    command {
        touch null.txt
        echo "one_line" > oneline.txt
        echo "" > emptyline.txt
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
    }
    output {
        Array[String] lines_null = read_lines("null.txt")
        Array[String] lines_oneline = read_lines("oneline.txt")
        Array[String] lines_emptyline = read_lines("emptyline.txt")
    }
}

workflow wf_test {
    call file_to_lines

    output {
        Int lines_null_count = length(file_to_lines.lines_null)
        Array[String] lines_oneline = file_to_lines.lines_oneline
        Array[String] lines_emptyline = file_to_lines.lines_emptyline
    }
}
