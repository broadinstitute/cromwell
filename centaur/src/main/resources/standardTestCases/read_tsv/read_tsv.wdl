task output_table {
  command {
    echo -e "col0\tcol1\tcol2"
    echo -e "a\tb\tc"
    echo -e "x\ty\tz"
  }
  output {
     Array[Array[String]] table = read_tsv(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task output_file_table {
  command {
    echo "first" > first
    echo "second" > second
    echo "third" > third
    echo "fourth" > fourth
    echo -e "first\tsecond"
    echo -e "third\tfourth"
  }
  output {
     Array[Array[File]] table = read_tsv(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task output_matrix {
  command {
    echo -e "0\t1\t2"
    echo -e "3\t4\t5"
    echo -e "6\t7\t8"
  }
  output {
     Array[Array[Int]] matrix = read_tsv(stdout())
  }
  runtime {
   docker: "ubuntu:latest"
  }
}

workflow test {
  call output_table
  call output_file_table
  call output_matrix
}
