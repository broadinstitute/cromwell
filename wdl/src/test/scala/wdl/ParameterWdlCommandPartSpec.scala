package wdl

import common.validation.Validation._
import wdl.command.ParameterCommandPart
import wdl.expression.NoFunctions
import wdl4s.parser.WdlParser.SyntaxError
import wom.types._
import wom.values._

import scala.util.Failure

class ParameterWdlCommandPartSpec extends WdlTest {
  val commandParameterWdl = "src/test/cases/command_parameters/test.wdl"

  "Command parameter WDL" should {
    val namespace = loadWdl("command_parameters/test.wdl")
    val task = namespace.tasks.find(_.name == "param_test") getOrElse {
      fail("task 'param_test' not found")
    }
    val paramsByName = task.commandTemplate.collect({ case p: ParameterCommandPart => p })

    s"Stringify the $${...} tags correctly" in {
      paramsByName.size shouldEqual 6
      paramsByName.head.toString shouldEqual s"$${a}"
      paramsByName(1).toString shouldEqual s"""$${"-p " + b}"""
      paramsByName(2).toString shouldEqual s"""$${sep="," c}"""
      paramsByName(3).toString shouldEqual s"""$${default="9" d}"""
      paramsByName(4).toString shouldEqual s"""$${sep="\\t" e}"""
      paramsByName(5).toString shouldEqual s"""$${true="--true" false="--false" f}"""
    }

    "raise exception if 'true' attribute is specified but 'false' is not" in {
      WdlNamespace.loadUsingSource(
        s"""task param_test {
            |  Boolean f
            |
            |  command <<<
            |  ./binary $${true="--true" f}
            |  >>>
            |}
            |
            |workflow wf {call param_test}
          """.stripMargin, None, None) match {
        case Failure(_: SyntaxError) => // expected
        case x => fail(s"Expecting a syntax error but got $x")
      }
    }

    "honor defined optional values during evaluation" in {
      val wdl =
        s"""
           |task t {
           |  String? some
           |  command {
           |     echo $${"hello " + some}
           |  }
           |}
        """.stripMargin

      val ns = WdlNamespace.loadUsingSource(wdl, None, None).get
      val task: WdlTask = ns.findTask("t").get
      val command = task.instantiateCommand(Map(task.declarations.head -> WomString("world")), NoFunctions).toTry.get.head
      command.commandString shouldBe "echo hello world"
    }
    
    "replace undefined values by their default value after evaluation" in {
      val wdl =
        s"""
          |task t {
          |  String? none
          |  command {
          |     echo $${"hello" + none}
          |  }
          |}
        """.stripMargin
      
      val ns = WdlNamespace.loadUsingSource(wdl, None, None).get
      val task = ns.findTask("t").get
      val command = task.instantiateCommand(Map(task.declarations.head -> WomOptionalValue.none(WomStringType)), NoFunctions).toTry.get.head
      command.commandString shouldBe "echo"
    }
  }
}
