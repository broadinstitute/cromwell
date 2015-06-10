package cromwell.util

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
}
