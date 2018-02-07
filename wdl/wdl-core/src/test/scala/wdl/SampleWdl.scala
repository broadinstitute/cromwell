package wdl

import java.nio.file.{Files, Path}

import wom.core.WorkflowSource
import wom.values.WomString

// FIXME: Figure out if anything can be removed from cromwell completely or pulled from here

trait SampleWdl extends TestFileUtil {
  def workflowSource(runtime: String = ""): WorkflowSource

  def createFileArray(base: Path): Unit = {
    createFile("f1", base, "line1\nline2\n")
    createFile("f2", base, "line3\nline4\n")
    createFile("f3", base, "line5\n")
    ()
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
    override def workflowSource(runtime: String = "") = sourceString()
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
  }

  object ThreeStep extends ThreeStepTemplate


  object NestedScatterWdl extends SampleWdl {
    override def workflowSource(runtime: String = "") =
      s"""
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
        |    python -c "print(len('$${B_in}'))"
        |  }
        |  output {
        |    Int B_out = read_int(stdout())
        |  }
        |}
        |
        |task C {
        |  Int C_in
        |  command {
        |    python -c "print($${C_in}*100)"
        |  }
        |  output {
        |    Int C_out = read_int(stdout())
        |  }
        |}
        |
        |task D {
        |  Array[Int] D_in
        |  command {
        |    python -c "print($${sep = '+' D_in})"
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
        |  scatter (item in A.A_out) { # scatter 0
        |    call B {input: B_in = item}
        |    call C {input: C_in = B.B_out}
        |    call E
        |    scatter (itemB in B.B_out) { # scatter 1
        |      call E as G
        |    }
        |    scatter (itemB in B.B_out) { # scatter 2
        |      call E as H
        |    }
        |  }
        |  scatter (item in A.A_out) { # scatter 3
        |    call E as F
        |  }
        |  call D {input: D_in = B.B_out}
        |}
      """.stripMargin
  }


  class ScatterWdl extends SampleWdl {
    val tasks = s"""task A {
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
      |    python -c "print(len('$${B_in}'))"
      |  }
      |  output {
      |    Int B_out = read_int(stdout())
      |  }
      |}
      |
      |task C {
      |  Int C_in
      |  command {
      |    python -c "print($${C_in}*100)"
      |  }
      |  output {
      |    Int C_out = read_int(stdout())
      |  }
      |}
      |
      |task D {
      |  Array[Int] D_in
      |  command {
      |    python -c "print($${sep = '+' D_in})"
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

    override def workflowSource(runtime: String = "") =
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
  }

  class DeclarationsWdl extends SampleWdl {
    override def workflowSource(runtime: String = "") =
      """task a {
        |  String foo = "notfoo"
        |  String bar = "bar"
        |  String foobar = foo + bar
        |  command {ps}
        |}
        |
        |task b {
        |  Int foo = 10
        |  Int bar = 2
        |  Int foobar = foo + bar
        |  Int foobar2 = foo + 2
        |  command {ps}
        |}
        |
        |task c {
        |  Int? foo = 3
        |  Array[Int]+ bar = [ 1, 2, 3 ]
        |  command {ps}
        |}
        |
        |workflow w {
        |  String foo = "foo"
        |  String bar = "bar"
        |  String foobar = foo + bar
        |  String? baz = "baz"
        |
        |  call a
        |  call a as a_prime
        |  call b
        |  call c
        |}
      """.stripMargin
  }

  object TaskDeclarationsWdl extends SampleWdl {
    override def workflowSource(runtime: String = "") =
      """
      |task t {
      |   String s
      |   command {
      |     echo ${s}
      |   }
      |   output {
      |     String o = s
      |     Array[Int] outputArray = [0, 1, 2]
      |   }
      |}
      |
      |task u {
      |   String a
      |   String b
      |   String c
      |   Int d
      |   String e = "e"
      |   String f
      |   String? g
      |   String? h
      |   String i
      |   File j
      |   Array[File] k
      |   String? l
      |
      |   command {
      |     echo ${a}
      |     echo ${b}
      |     echo ${c}
      |     echo ${d}
      |     echo ${e}
      |     echo ${f}
      |     echo ${g}
      |     echo ${h}
      |     echo ${i}
      |   }
      |}
      |
      |workflow wf {
      | String workflowDeclarationFromInput
      | String workflowDeclaration = "b"
      | Array[File] files = ["a", "b", "c"]
      |
      | call t as t2 {input: s = "hey" }
      |
      | scatter (i in t2.outputArray) {
      |   call t {input: s = "c"}
      |   if (true) {
      |     call t as t3 {input: s = "c"}
      |   }
      |   call u as v {input: a = workflowDeclarationFromInput,
      |                       b = workflowDeclaration,
      |                       c = t.o,
      |                       d = i,
      |                       i = "${workflowDeclaration}",
      |                       k = files,
      |                       l = t3.o }
      | }
      |}
    """.
        stripMargin

    val workflowInputs = Map(
      "wf.workflowDeclarationFromInput" -> WomString("a"),
      "wf.v.f" -> WomString("f"),
      "wf.v.g" -> WomString("g"),
      "wf.v.j" -> WomString("j")
    )
  }
}
