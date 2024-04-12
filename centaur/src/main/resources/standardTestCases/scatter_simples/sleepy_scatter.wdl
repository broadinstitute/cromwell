version 1.0

task genSeq {
  input {
    Int count
  }
  command <<<
    echo "~{count}"
    seq 1 ~{count}
  >>>

  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }

  output {
    Array[Int] seq = read_lines(stdout())
  }
}

task nap {
  input {
    Int secs
    Int napIdx
  }
  command <<<
    echo "Nap ${napIdx} commencing..."
    sleep ${secs}
    echo "Ah, that was refreshing!"
  >>>

  runtime {
    docker: "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
  }

  output {
    String napReport = stdout()
  }
}

workflow w {
  input {
    Int shardCount = 3
    Int sleepyTime = 5
  }

  call genSeq {
    input: count = shardCount
  }
  scatter (idx in genSeq.seq) {
    call nap {
      input:
        secs = sleepyTime,
        napIdx = idx
    }
  }
  output {
    String done = "done!"
  }
}
