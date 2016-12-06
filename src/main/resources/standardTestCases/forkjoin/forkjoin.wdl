##
# Checks a simple branch and join operation.
# We start with a task, branch into two parallel executions, and then rejoin to calculate the result.
##

task mkFile {
  command {
    for i in `seq 1 1000`
    do
      echo $i
    done
    sleep 2
  }
  output {
    File numbers = stdout()
  }
  runtime {docker: "ubuntu:latest"}
}

task grep {
  String pattern
  File in_file
  command {
    grep '${pattern}' ${in_file} | wc -l
    sleep 2
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task wc {
  File in_file
  command {
    cat ${in_file} | wc -l
    sleep 2
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task join {
  Int grepCount
  Int wcCount
  command {
    expr ${wcCount} / ${grepCount}
    sleep 2
  }
  output {
    Int proportion = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

workflow forkjoin {
  call mkFile
  call grep { input: in_file = mkFile.numbers }
  call wc { input: in_file=mkFile.numbers }
  call join { input: wcCount = wc.count, grepCount = grep.count }
  output {
    join.proportion
  }
}
