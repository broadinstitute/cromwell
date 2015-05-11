package cromwell

import cromwell.binding.WdlBinding

trait WdlThreeStepFixture {
  val binding = WdlBinding.process(
    """
      |task ps {
      |  command {
      |    ps
      |  }
      |  output {
      |    File procs = "stdout"
      |  }
      |}
      |
      |task cgrep {
      |  command {
      |    grep '${pattern}' ${File in_file} | wc -l
      |  }
      |  output {
      |    Int count = read_int("stdout")
      |  }
      |}
      |
      |task wc {
      |  command {
      |    wc -l ${File in_file}
      |  }
      |  output {
      |    Int count = read_int("stdout")
      |  }
      |}
      |
      |workflow three_step {
      |  call ps
      |  call cgrep {
      |    input: in_file=ps.procs
      |  }
      |  call wc {
      |    input: in_file=ps.procs
      |  }
      |}
      |
    """.stripMargin)
}
