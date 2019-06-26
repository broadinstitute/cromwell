version 1.0

workflow cached_inputs {
    Array[Int] one_to_ten = [1,2,3,4,5,6,7,8,9,10]

    call ten_lines

    scatter (x in one_to_ten) {
        call read_line {
            input:
                file=ten_lines.text,
                line_number=x
        }
    }
    output {
        Array[String] lines = read_line.line
    }
}

task ten_lines {
    command {
        echo "Line 1
        Line 2
        Line 3
        Line 4
        Line 5
        Line 6
        Line 7
        Line 8
        Line 9
        Line 10" > outfile.txt
    }
    output {
        File text = "outfile.txt"
    }
    runtime {
      docker: "ubuntu:latest"
    }
}

task read_line {
    input {
        File file
        Int line_number
    }
    command {
        sed -n ~{line_number}p ~{file}
    }
    output {
        String line = read_string(stdout())
    }
    runtime {
      docker: "ubuntu:latest"
    }
}