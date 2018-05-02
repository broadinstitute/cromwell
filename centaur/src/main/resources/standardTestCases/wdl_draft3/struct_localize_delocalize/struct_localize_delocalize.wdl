version 1.0

struct FileContainer {
  File bam
  File index
}

workflow struct_localize_delocalize {
  call struct_delocalize
  call struct_localize { input: container = struct_delocalize.fc }

  scatter(s in struct_localize.found_files) {
    String b = basename(s)
  }

  output {
    Array[String] result = b
  }
}

task struct_delocalize {
  command <<<
    echo "BAM BAM BAM BAM BAM" > test.bam
    echo "BAI BAI BAI BAI BAI" > test.bam.bai
  >>>
  output {
    FileContainer fc = object {
      bam: "test.bam",
      index: "test.bam.bai"
    }
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

task struct_localize {
  input {
    FileContainer container
  }
  command <<<
    # Should find the test.bam and test.bam.bai
    ls ~{container.bam}* > bam_files
  >>>
  output {
    Array[String] found_files = read_lines("bam_files")
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
