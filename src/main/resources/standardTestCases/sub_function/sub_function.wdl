task make_a_file {
  command {
    echo "content" > myfilename.bam
  }
  output {
    File a_file = "myfilename.bam"
  }
  runtime {docker: "ubuntu:latest"}
}

task sub {
    String myBamString = "myfilename.bam"
    File myBamFile
    String swappedStr = sub(myBamString, ".bam$", ".txt")
    # At this point myBamFile is not localized so the protocol and path must be removed.
    String swappedFile = sub(sub(sub(myBamFile, "^[a-z0-9]+:", ""), "/.*/",""), ".bam$", ".txt")

    command {
      echo ${sub(myBamString, ".bam$", ".txt")}
      echo "content" > ${swappedFile}
    }

    runtime {
        docker: "ubuntu:latest"
    }

    output {
      String stdout = read_string(stdout())
      String swappedContent = read_string(swappedFile)
      String o2 = swappedStr
    }
}

workflow wf {
    call make_a_file
    call sub {
        input: myBamFile = make_a_file.a_file
    }
}

