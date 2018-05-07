version 1.0

##
# Checks a simple branch and join operation.
# We start with a task, branch into two parallel executions, and then rejoin to calculate the result.
##

workflow forkjoin {
  call mkFile

  call grep { input: in_file = mkFile.numbers }
  call wc { input: in_file=mkFile.numbers }

  call join { input: wcCount = wc.count, grepCount = grep.count }

  output {
    Int joined_proportion = join.proportion
  }
}

task mkFile {
  command <<<
    for i in `seq 1 1000`
    do
      echo $i
    done
  >>>
  output {
    File numbers = stdout()
  }
  runtime {docker: "ubuntu:latest"}
}

task grep {
  input {
    String pattern
    File in_file
  }
  command <<<
    [ -f "~{in_file}" ] && grep '~{pattern}' ~{in_file} | wc -l
  >>>
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task wc {
  input {
    File in_file
  }
  command <<<
    [ -f "~{in_file}" ] && cat ~{in_file} | wc -l
  >>>
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task join {
  input {
    Int grepCount
    Int wcCount
  }
  command <<<
    expr ~{wcCount} / ~{grepCount}
  >>>
  output {
    Int proportion = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}
