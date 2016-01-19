package wdl4s

import wdl4s.parser.WdlParser.SyntaxError
import org.scalatest.{FlatSpec, Matchers}

class SyntaxErrorSpec extends FlatSpec with Matchers {
  val psTaskWdl = """
      |task ps {
      |  command {
      |    ps
      |  }
      |  output {
      |    File procs = stdout()
      |  }
      |}""".stripMargin

  val cgrepTaskWdl = """
     |task cgrep {
     |  String pattern
     |  File in_file
     |  command {
     |    grep '${pattern}' ${in_file} | wc -l
     |  }
     |  output {
     |    Int count = read_int(stdout())
     |  }
     |}""".stripMargin

  def resolver(importUri: String): WdlSource = {
    importUri match {
      case "ps" => psTaskWdl
      case "cgrep" => cgrepTaskWdl
      case _ => throw new RuntimeException(s"Can't resolve $importUri")
    }
  }

  private def expectError(wdl: String) = {
    try {
      val namespace = WdlNamespace.load(wdl, resolver _)
      fail("Exception expected")
    } catch {
      case x: SyntaxError => // expected
    }
  }

  "WDL syntax checker" should "detect error when call references bad input" in {
    expectError("""
      |task ps {
      |  command {
      |    ps
      |  }
      |  output {
      |    File procs = stdout()
      |  }
      |}
      |task cgrep {
      |  String pattern
      |  File in_file
      |  command {
      |    grep '${pattern}' ${in_file} | wc -l
      |  }
      |  output {
      |    Int count = read_int(stdout())
      |  }
      |}
      |workflow three_step {
      |  call ps
      |  call cgrep {
      |    input: BADin_file=ps.procs
      |  }
      |}""".stripMargin)
  }

  it should "detect error when call references bad task" in {
    expectError("""
      |task ps {
      |  command {
      |    ps
      |  }
      |  output {
      |    File procs = stdout()
      |  }
      |}
      |task cgrep {
      |  String pattern
      |  File in_file
      |  command {
      |    grep '${pattern}' ${in_file} | wc -l
      |  }
      |  output {
      |    Int count = read_int(stdout())
      |  }
      |}
      |workflow three_step {
      |  call ps
      |  call cgrepBAD {
      |    input: in_file=ps.procs
      |  }
      |}""".stripMargin)
  }
  it should "detect error when more than one workflow is defined" in {
    expectError("""
        |task ps {
        |  command {
        |    ps
        |  }
        |  output {
        |    File procs = stdout()
        |  }
        |}
        |task cgrep {
        |  File in_file
        |  String pattern
        |  command {
        |    grep '${pattern}' ${in_file} | wc -l
        |  }
        |  output {
        |    Int count = read_int(stdout())
        |  }
        |}
        |workflow three_step {
        |  call ps
        |  call cgrep {
        |    input: in_file=ps.procs
        |  }
        |}
        |workflow BAD {}""".stripMargin)
  }
  it should "detect error when namespace and task have the same name" in {
    expectError("""
        |import "ps" as ps
        |task ps {command {ps}}
        |workflow three_step {
        |  call ps
        |}
        |""".stripMargin)
  }
  it should "detect error when namespace and workflow have the same name" in {
    expectError("""
        |import "ps" as ps
        |workflow ps {
        |  call ps
        |}
        |""".stripMargin)
  }
  it should "detect error when two tasks have the same name" in {
    expectError("""
        |import "ps"
        |task ps {command {ps}}
        |workflow three_step {
        |  call ps
        |}
        |""".stripMargin)
  }
  it should "detect error a MemberAccess references a non-existent input on a task" in {
    expectError("""
        |import "ps"
        |import "cgrep"
        |workflow three_step {
        |  call ps
        |  call cgrep {
        |    input: pattern=ps.BAD
        |  }
        |}
        |""".stripMargin)
  }
  it should "detect error a MemberAccess references a non-existent left-hand side" in {
    expectError("""
        |import "ps"
        |import "cgrep"
        |workflow three_step {
        |  call ps
        |  call cgrep {
        |    input: pattern=psBAD.procs
        |  }
        |}
        |""".stripMargin)
  }
  it should "detect unexpected EOF" in {
    expectError("workflow")
  }
  it should "detect unexpected symbol" in {
    expectError("workflow foo workflow")
  }
  it should "detect extraneous symbols" in {
    expectError("workflow foo {}}")
  }
  it should "detect when two call definitions have the same name" in {
    expectError(
      """
        |task x {
        |  command { ps }
        |}
        |
        |workflow wf {
        |  call x
        |  call x
        |}
      """.stripMargin)
  }
  it should "detect when there are more than one 'input' sections in a call" in {
    expectError(
      """
        |task x {
        |  String a
        |  String b
        |  command {  ./script ${a} ${b} }
        |}
        |
        |workflow wf {
        |  call x {
        |    input: a = "a"
        |    input: b = "b"
        |  }
        |}
      """.stripMargin)
  }
  it should "detect when there are two workflows defined in a WDL file" in {
    expectError(
      """workflow w {}
        |workflow x {}
      """.stripMargin)
  }
  it should "detect when a Map does not have two parameterized types" in {
    expectError(
      """workflow w {
        |  Map[Int] i
        |}
      """.stripMargin)
  }
  it should "detect when task output section declares an output with incompatible types" in {
    expectError(
      """task a {
        |  command { ./script }
        |  output {
        |    Array[String] x = "bad value"
        |  }
        |}
        |
        |workflow w {
        |  call a
        |}
      """.stripMargin)
  }
  it should "detect when task output section declares an output with incompatible types 2" in {
    expectError(
      """task a {
        |  command { ./script }
        |  output {
        |    Int x = "bad value"
        |  }
        |}
        |
        |workflow w {
        |  call a
        |}
      """.stripMargin)
  }
  it should "detect when a meta or parameter_meta value is not a string" in {
    expectError(
      """task a {
        |  command { ./script }
        |  meta {
        |    foo: 1+1
        |  }
        |}
        |
        |workflow w {
        |  call a
        |}
      """.stripMargin)
  }
  it should "detect when a two command sections are specified" in {
    expectError(
      """task a {
        |  command { ./script }
        |  command { ps }
        |}
        |
        |workflow w {
        |  call a
        |}
      """.stripMargin)
  }
  it should "detect when expressions within a command block don't reference actual variables" in {
    expectError(
      """task a {
        |  Int x
        |  command { ./script ${x+y} }
        |}
        |
        |workflow w {
        |  call a
        |}
      """.stripMargin)
  }
  it should "detect when a task declaration references a variable that doesn't exist" in {
    expectError(
      """task a {
        |  Int x
        |  Int y = x + z
        |  command { ./script ${x} }
        |}
        |
        |workflow w {
        |  call a {input: x=5}
        |}
      """.stripMargin)
  }
  it should "detect when there's more than one command in a task" in {
    expectError(
      """task a {
        |  Int x
        |  command { ./script ${x} }
        |  command { ps }
        |}
        |
        |workflow w {
        |  call a {input: x=5}
        |}
      """.stripMargin)
  }
  it should "detect duplicate variable declarations (1)" in {
    expectError(
      """task inputOops {
        |  Int a = 5
        |  Int a = 10
        |  command { echo ${a} }
        |}
        |workflow a { call inputOops }
      """.stripMargin)
  }
  it should "detect duplicate variable declarations (2)" in {
    expectError(
      """task outputOops {
        |  command { echo 5 }
        |  output {
        |    Int b = 5
        |    Int b = 10
        |  }
        |}
        |workflow b { call outputOops }
      """.stripMargin)
  }
  it should "detect duplicate variable declarations (3)" in {
    expectError(
      """task inputOutputOops {
        |  Int c = 5
        |  command { echo 5 }
        |  output {
        |    Int c = 10
        |  }
        |}
        |workflow c { call inputOutputOops }
      """.stripMargin)
  }
}

