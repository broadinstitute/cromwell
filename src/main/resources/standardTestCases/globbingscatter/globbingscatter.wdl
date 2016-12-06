task globber {
  Int count
  command <<<
    for i in `seq 1 ${count}`
    do
      mkdir out-$i
      echo "globbing is my number $i best hobby" > out-$i/$i.txt
    done
    sleep 2
  >>>
  output {
    Array[File] out = glob("out-*/*.txt")
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

task catter {
  File in_file
  command <<<
    cat ${in_file}
    sleep 2
  >>>
  output {
    File result = stdout()
  }
  runtime {docker:"ubuntu:latest"}
}

workflow globbingscatter {
  call globber { input: count=5 }
  call combiner as combiner1 { input: in_file=globber.out }
  scatter (x in globber.out) {
    call catter as catter1 { input: in_file=x }
    call catter as catter2 { input: in_file=catter1.result }
  }
  call combiner as combiner2 { input: in_file=catter2.result }
  output {
     combiner2.result
  }
}
