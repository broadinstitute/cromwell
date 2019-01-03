version 1.0
#
# Goal here is to split up the input file into files of 1 line each (in the prepare) then in parallel call wc -w on each newly created file and count the words into another file then in the gather, su
# the word-count for the file
#
# splits each line into a file with the name temp_?? (shuffle)
task do_prepare {
    Array[String] lines = [
      "the",
      "total number",
      "of words in this",
      "text file is 11"
    ]
    command {
        split -l 1 ~{write_lines(lines)} temp_
    }
    output {
        Array[File] split_files = ["temp_aa", "temp_ab", "temp_ac", "temp_ad"]
    }
    runtime {
      docker: "ubuntu:latest"
    }
}
# count the number of words in the input file, writing the count to an output file overkill in this case, but simulates a real scatter-gather that would just return an Int (map)
task do_scatter {
    input {
      String salt
      File input_file
    }
    command {
        # ~{salt}
        wc -w ~{input_file} > output.txt
    }
    output {
        File count_file = "output.txt"
    }
    runtime {
      docker: "ubuntu:latest"
    }
}
# aggregate the results back together (reduce)
task do_gather {
    input {
      Array[File] input_files
    }
    command <<<
        cat ~{sep = ' ' input_files} | awk '{s+=$1} END {print s}'
    >>>
    output {
      Int sum = read_int(stdout())
    }
    runtime {
      docker: "ubuntu:latest"
    }
}

workflow sc_test {
    call do_prepare
    scatter(f in do_prepare.split_files) {
      call do_scatter {
            input: input_file = f
        }
    }
    call do_gather {
        input: input_files = do_scatter.count_file
    }
}
