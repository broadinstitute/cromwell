package cromwell.util

import java.io.{File, FileWriter}

import cromwell.binding._

trait SampleWdl {
  def wdlSource(runtime: String = ""): WdlSource

  val rawInputs: WorkflowRawInputs

  def wdlJson: WdlJson = {
    "{" + rawInputs.collect { case (k, v) => s""" "$k": "$v"""" }.mkString(",\n") + "}"
  }
}

object SampleWdl {

  object HelloWorld extends SampleWdl {
    override def wdlSource(runtime: String = "") =
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
    val rawInputs = Map(Addressee -> "world")
    val OutputKey = "hello.hello.salutation"
    val OutputValue = "Hello world!\n"
  }

  object GoodbyeWorld extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        |task goodbye {
        |  command {
        |    sh -c "exit 1"
        |  }
        |  output {
        |    String out = read_string(stdout())
        |  }
        |}
        |
        |workflow goodbye {
        |  call goodbye
        |}
      """.stripMargin

    val rawInputs = Map.empty[String, Any]
    val OutputKey = "goodbye.goodbye.out"
  }

  object Incr extends SampleWdl {
    override def wdlSource(runtime: String = "") =
    """
      |task incr {
      |  command {
      |    echo $((${Int val} + 1))
      |  }
      |  output {
      |    Int out = read_int(stdout())
      |  }
      |}
      |
      |workflow incr {
      |  call incr
      |}
    """.stripMargin

    override val rawInputs: WorkflowRawInputs = Map("incr.incr.val" -> "1")
  }

  object SubtractionWorkflow {
    val WdlSource =
      """
        |task a {
        |  command { echo '${message}' }
        |  output {
        |    String message = read_string(stdout())
        |    Int constant = 100
        |  }
        |}
        |task b {
        |  command { echo '${message} - ${Int integer}' }
        |}
        |workflow wf {
        |  call a
        |  call b {
        |    input: message=a.message, integer=a.constant - 75
        |  }
        |}
      """.stripMargin
  }

  object CoercionNotDefined extends SampleWdl {
    override def wdlSource(runtime: String = "") = {
      """
        |task summary {
        |  command {
        |    ~/plink --bfile ${bfile} --missing --hardy --out foo --allow-no-sex
        |  }
        |  output {
        |    File hwe = "foo.hwe"
        |    File log = "foo.log"
        |    File imiss = "foo.imiss"
        |    File lmiss = "foo.lmiss"
        |  }
        |  meta {
        |    author: "Jackie Goldstein"
        |    email: "jigold@broadinstitute.org"
        |  }
        |}
        |
        |
        |workflow test1 {
        |  String bfile
        |  call summary {
        |     input: bfile=bfile
        |  }
        |}
      """.stripMargin
    }

    override val rawInputs: WorkflowRawInputs = Map("test1.bfile" -> "data/example1")
  }

  object ThreeStep extends SampleWdl {
    override def wdlSource(runtime: String = "") =
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

    val PatternKey = "three_step.cgrep.pattern"
    override val rawInputs = Map(PatternKey -> "...")
  }

  object ThreeStepLargeJson extends SampleWdl {
    override def wdlSource(runtime: String = "") = ThreeStep.wdlSource(runtime)
    override lazy val rawInputs = Map(ThreeStep.PatternKey -> "." * 10000)
  }

  object CannedThreeStep extends SampleWdl {
    val CannedProcessOutput =
      """
        |USER              PID  %CPU %MEM      VSZ    RSS   TT  STAT STARTED      TIME COMMAND
        |joeblaux         10344   4.5 20.6  7011548 3454616  ??  S    Mon06AM 275:26.10 /Applications/IntelliJ IDEA 14.app/Contents/MacOS/idea
        |joeblaux           883   2.2  0.5  2716336  85768   ??  S    Sun08AM   0:52.64 /Applications/Utilities/Terminal.app/Contents/MacOS/Terminal
        |_coreaudiod        577   1.9  0.1  2522428   9572   ??  Ss   Sun08AM  55:39.69 /usr/sbin/coreaudiod
        |_windowserver      130   1.6  1.9  4535136 319588   ??  Ss   Sun08AM 148:39.39 /System/Library/Frameworks/ApplicationServices.framework/Frameworks/CoreGraphics.framework/Resources/WindowServer -daemon
        |joeblaux          6440   1.4  2.2  4496164 362136   ??  S    Sun09PM  74:29.40 /Applications/Google Chrome.app/Contents/MacOS/Google Chrome
      """.stripMargin.trim

    override def wdlSource(runtime: String): WdlSource =
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

    private def createCannedPsFile: File = {
      val file = File.createTempFile("canned_ps", ".out")
      val writer = new FileWriter(file)
      writer.write(CannedProcessOutput)
      writer.flush()
      writer.close()
      file
    }

    override val rawInputs = {
      import InputKeys._
      Map(
        Pattern -> "joeblaux",
        DummyPsFile -> createCannedPsFile.getAbsolutePath,
        DummyPs2File -> createCannedPsFile.getAbsolutePath,
        DummyPs3File -> createCannedPsFile.getAbsolutePath)
    }
  }

  object CannedFilePassing extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        |task ps {
        |  command {
        |    (echo "x"; echo "y"; echo "z") > myfile.txt
        |  }
        |  output {
        |    File procs = "myfile.txt"
        |  }
        |  runtime {
        |    docker: "ubuntu:latest"
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
        |  runtime {
        |    docker: "ubuntu:latest"
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
        |  runtime {
        |    docker: "ubuntu:latest"
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

    override val rawInputs = Map(ThreeStep.PatternKey -> "x")

    val Expectations = Map("three_step.wc.count" -> 3, "three_step.cgrep.count" -> 1)
  }

  object InputKeys {
    val Pattern = "three_step.cgrep.pattern"
    // ps2 and ps3 are not part of the core three-step workflow, but are intended to flush out issues
    // with incorrectly starting multiple copies of cgrep and wc calls due to race conditions.
    val DummyPsFile = "three_step.ps.dummy_ps_file"
    val DummyPs2File = "three_step.ps2.dummy_ps_file"
    val DummyPs3File = "three_step.ps3.dummy_ps_file"
  }

  object CurrentDirectory extends SampleWdl {
    override def wdlSource(runtime: String): String =
      """
        |task whereami {
        |  command {
        |    pwd
        |  }
        |  output {
        |    String pwd = read_string(stdout())
        |  }
        |}
        |
        |workflow whereami {
        |  call whereami
        |}
      """.stripMargin

    override val rawInputs: Map[String, Any] = Map.empty
  }

}
