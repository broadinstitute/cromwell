package cromwell.util

import cromwell.binding._

object SampleWdl {
  object HelloWorld {
    val WdlSource =
      """
        |task hello {
        |  command {
        |    echo "Hello ${addressee}!"
        |  }
        |  output {
        |    String salutation = read_string(stdout())
        |  }
        |}
        |
        |workflow hello {
        |  call hello
        |}
      """.stripMargin

    val Addressee = "hello.hello.addressee"
    val RawInputs =  Map(Addressee -> "world")
    val OutputKey = "hello.hello.salutation"
    val OutputValue = "Hello world!\n"
    val JsonInputs = s""" { "$Addressee" : "world" } """
  }

  object ThreeStep {
      val WdlSource =
        """
        |task ps {
        |  command {
        |    ps
        |  }
        |  output {
        |    File procs = stdout()
        |  }
        |}
        |
        |task cgrep {
        |  command {
        |    grep '${pattern}' ${File in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |}
        |
        |task wc {
        |  command {
        |    cat ${File in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
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
      """.stripMargin

    val RawInputs =  Map("three_step.cgrep.pattern" -> "...")
  }

  object FauxThreeStep {
    def wdlSource(runtime: String): WdlSource =
      """
        |task ps {
        |  command {
        |    cat ${File dummy_ps_file}
        |  }
        |  output {
        |    File procs = stdout()
        |  }
        |  RUNTIME
        |}
        |
        |task cgrep {
        |  command {
        |    grep '${pattern}' ${File in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |  RUNTIME
        |}
        |
        |task wc {
        |  command {
        |    cat ${File in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |  RUNTIME
        |}
        |
        |workflow three_step {
        |  call ps
        |  call ps as ps2
        |  call ps as ps3
        |  call cgrep {
        |    input: in_file=ps.procs
        |  }
        |  call wc {
        |    input: in_file=ps.procs
        |  }
        |}
      """.stripMargin.replaceAll("RUNTIME", runtime)


    object InputKeys {
      val Pattern = "three_step.cgrep.pattern"
      // ps2 and ps3 are not part of the core three-step workflow, but are intended to flush out issues
      // with incorrectly starting multiple copies of cgrep and wc calls due to race conditions.
      val DummyPsFile = "three_step.ps.dummy_ps_file"
      val DummyPs2File = "three_step.ps2.dummy_ps_file"
      val DummyPs3File = "three_step.ps3.dummy_ps_file"
    }
  }
}
