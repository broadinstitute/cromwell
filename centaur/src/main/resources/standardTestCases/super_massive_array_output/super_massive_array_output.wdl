task StringSpam {
  String sample_name = "DeliciousStringSpam"
  Int index
  Int spamWidth

  command <<<
    mkdir ${sample_name}_${index}
    for i in `seq 0 ${spamWidth}`
    do
      echo "Yum yum yum what a delicious piece of spam. This is long to increase the final output size. And a little more text just to increase that size just a little bit more... sample_name_${index}_$i" >> file.txt
    done
  >>>
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Array[String] outArray = read_lines("file.txt")
  }
}

task MatrixRotation {
  Array[Array[String]] input_matrix

    command <<<
      python <<CODE
      import csv
      import sys
      with open('${write_tsv(input_matrix)}') as tsv_in:
        input_matrix = [line.strip().split('\t') for line in tsv_in]
        final_matrix = [["" for x in range(len(input_matrix))] for y in range(len(input_matrix[0]))]
        for x in range(len(input_matrix)):
          for y in range(len(input_matrix[0])):
            final_matrix[y][x] = input_matrix[x][y]
      with open("output.tsv", "w") as tsv_out:
        # lineterminator is a workaround for a cromwell bug that doesnt allow for '\r\n' line endings which this outputs by default
        writer = csv.writer(tsv_out, delimiter='\t', lineterminator='\n')
        writer.writerows(final_matrix)
      CODE
      >>>

    runtime {
      docker: "python:2.7"
      memory: "1 GB"
      preemptible: 3
    }

    output {
        File out = "output.tsv"
        Array[Array[String]] output_matrix = read_tsv(out)
    }
}

workflow DeliciousFileSpam {

  Array[Int] indexing_list = range(500)

  scatter (idx in indexing_list) {
    call StringSpam { input: index = idx, spamWidth = 900 }
  }
}