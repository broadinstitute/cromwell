version 1.0

task genSeq {
  input {
    Int count
  }
  command <<<
    echo "~{count}"
    seq 1 ~{count}
  >>>
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

  output {
    String napReport = stdout()
  }
}

workflow w {
  Int shardCount = 5
  Int sleepyTime = 3

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
