package wdl4s

import java.nio.file.{Files, Path}

import scala.language.postfixOps

// FIXME: Figure out if anything can be removed from cromwell completely or pulled from here

trait SampleWdl extends TestFileUtil {
  def wdlSource(runtime: String = ""): WdlSource
  def rawInputs: WorkflowRawInputs

  def createFileArray(base: Path): Unit = {
    createFile("f1", base, "line1\nline2\n")
    createFile("f2", base, "line3\nline4\n")
    createFile("f3", base, "line5\n")
  }

  def cleanupFileArray(base: Path) = {
    deleteFile(base.resolve("f1"))
    deleteFile(base.resolve("f2"))
    deleteFile(base.resolve("f3"))
  }

  def deleteFile(path: Path) = Files.delete(path)
}

object SampleWdl {
  trait ThreeStepTemplate extends SampleWdl {
    override def wdlSource(runtime: String = "") = sourceString()
    private val outputSectionPlaceholder = "OUTPUTSECTIONPLACEHOLDER"
    def sourceString(outputsSection: String = "") = {
      val withPlaceholders =
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
        |  String pattern
        |  File in_file
        |
        |  command {
        |    grep '${pattern}' ${in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |}
        |
        |task wc {
        |  File in_file
        |  command {
        |    cat ${in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |}
        |
        |workflow three_step {
        |  call ps
        |  call cgrep {
        |    input: in_file = ps.procs
        |  }
        |  call wc {
        |    input: in_file = ps.procs
        |  }
        |  """ + outputSectionPlaceholder + """
        |}
        |
        """
      withPlaceholders.stripMargin.replace(outputSectionPlaceholder, outputsSection)
    }

    val PatternKey ="three_step.cgrep.pattern"
    override lazy val rawInputs = Map(PatternKey -> "...")
  }

  object ThreeStep extends ThreeStepTemplate


  object NestedScatterWdl extends SampleWdl {
    override def wdlSource(runtime: String = "") =
      """
        task A {
        |  command {
        |    echo -n -e "jeff\nchris\nmiguel\nthibault\nkhalid\nscott"
        |  }
        |  output {
        |    Array[String] A_out = read_lines(stdout())
        |  }
        |}
        |
        |task B {
        |  String B_in
        |  command {
        |    python -c "print(len('${B_in}'))"
        |  }
        |  output {
        |    Int B_out = read_int(stdout())
        |  }
        |}
        |
        |task C {
        |  Int C_in
        |  command {
        |    python -c "print(${C_in}*100)"
        |  }
        |  output {
        |    Int C_out = read_int(stdout())
        |  }
        |}
        |
        |task D {
        |  Array[Int] D_in
        |  command {
        |    python -c "print(${sep = '+' D_in})"
        |  }
        |  output {
        |    Int D_out = read_int(stdout())
        |  }
        |}
        |
        |task E {
        |  command {
        |    python -c "import random; print(random.randint(1,100))"
        |  }
        |  output {
        |    Int E_out = read_int(stdout())
        |  }
        |}
        |
        |workflow w {
        |  call A
        |  scatter (item in A.A_out) {
        |    call B {input: B_in = item}
        |    call C {input: C_in = B.B_out}
        |    call E
        |    scatter (itemB in B.B_out) {
        |     call E as G
        |    }
        |    scatter (itemB in B.B_out) {
        |     call E as H
        |    }
        |  }
        |  scatter (item in A.A_out) {
        |    call E as F
        |  }
        |  call D {input: D_in = B.B_out}
        |}
      """.stripMargin
    override lazy val rawInputs = Map("" -> "...")
  }


  class ScatterWdl extends SampleWdl {
    val tasks = """task A {
      |  command {
      |    echo -n -e "jeff\nchris\nmiguel\nthibault\nkhalid\nscott"
      |  }
      |  output {
      |    Array[String] A_out = read_lines(stdout())
      |  }
      |}
      |
      |task B {
      |  String B_in
      |  command {
      |    python -c "print(len('${B_in}'))"
      |  }
      |  output {
      |    Int B_out = read_int(stdout())
      |  }
      |}
      |
      |task C {
      |  Int C_in
      |  command {
      |    python -c "print(${C_in}*100)"
      |  }
      |  output {
      |    Int C_out = read_int(stdout())
      |  }
      |}
      |
      |task D {
      |  Array[Int] D_in
      |  command {
      |    python -c "print(${sep = '+' D_in})"
      |  }
      |  output {
      |    Int D_out = read_int(stdout())
      |  }
      |}
      |
      |task E {
      |  command {
      |    python -c "print(9)"
      |  }
      |  output {
      |    Int E_out = read_int(stdout())
      |  }
      |}
    """.stripMargin

    override def wdlSource(runtime: String = "") =
      s"""$tasks
        |
        |workflow w {
        |  call A
        |  scatter (item in A.A_out) {
        |    call B {input: B_in = item}
        |    call C {input: C_in = B.B_out}
        |    call E
        |  }
        |  call D {input: D_in = B.B_out}
        |}
      """.stripMargin

    override lazy val rawInputs = Map.empty[String, String]
  }
}
