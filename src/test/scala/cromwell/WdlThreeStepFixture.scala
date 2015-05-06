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
          |    file procs = "stdout"
          |  }
          |}
          |
          |task cgrep {
          |  command {
          |    grep '${pattern}' ${file in_file} | wc -l
          |  }
          |  output {
          |    int count = read_int("stdout")
          |  }
          |}
          |
          |task wc {
          |  command {
          |    wc -l ${file in_file}
          |  }
          |  output {
          |    int count = read_int("stdout")
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
