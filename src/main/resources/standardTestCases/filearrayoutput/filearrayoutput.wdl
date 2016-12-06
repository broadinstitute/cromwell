task createFileArray {
  command <<<
    mkdir out
    echo "hullo" > out/hello.txt
    echo "buh-bye" > out/ciao.txt
    sleep 2
  >>>
  output {
    Array[File] out = [ "out/hello.txt", "out/ciao.txt" ]
  }
  runtime {docker:"ubuntu:latest"}
}

task combiner {
  Array[File] in_file
  command <<<
    cat ${sep=' ' in_file}
    sleep 2
  >>>
  output {
    String result = read_string(stdout())
  }
  runtime {docker:"ubuntu:latest"}
}

workflow filearrayoutput {
    call createFileArray
    call combiner { input: in_file = createFileArray.out }
    output {
       combiner.result
    }
}
