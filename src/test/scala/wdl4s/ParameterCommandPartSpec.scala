package wdl4s

import better.files._
import wdl4s.command.ParameterCommandPart
import wdl4s.expression.{NoFunctions, PureStandardLibraryFunctions}
import wdl4s.parser.WdlParser.SyntaxError
import wdl4s.types.WdlStringType
import wdl4s.values.{WdlOptionalValue, WdlString}

import scala.util.{Failure, Success, Try}

class ParameterCommandPartSpec extends WdlTest {
  val commandParameterWdl = "src/test/cases/command_parameters/test.wdl"

  commandParameterWdl should {
    val namespace = loadWdlFile(File(commandParameterWdl))
    val task = namespace.tasks.find(_.name == "param_test") getOrElse {
      fail("task 'param_test' not found")
    }
    val paramsByName = task.commandTemplate.collect({ case p: ParameterCommandPart => p })

    "Stringify the ${...} tags correctly" in {
      paramsByName.size shouldEqual 6
      paramsByName.head.toString shouldEqual "${a}"
      paramsByName(1).toString shouldEqual "${\"-p \" + b}"
      paramsByName(2).toString shouldEqual "${sep=\",\" c}"
      paramsByName(3).toString shouldEqual "${default=\"9\" d}"
      paramsByName(4).toString shouldEqual "${sep=\"\\t\" e}"
      paramsByName(5).toString shouldEqual "${true=\"--true\" false=\"--false\" f}"
    }

    "raise exception if 'true' attribute is specified but 'false' is not" in {
      Try(
        WdlNamespace.loadUsingSource(
          """task param_test {
            |  Boolean f
            |
            |  command <<<
            |  ./binary ${true="--true" f}
            |  >>>
            |}
            |
            |workflow wf {call param_test}
          """.stripMargin, None, None)
      ) match {
        case Failure(s: SyntaxError) => // expected
        case _ => fail("Expecting a syntax error")
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

      val ns = WdlNamespace.loadUsingSource(wdl, None, None)
      val task: Task = ns.findTask("t").get
      val command = task.instantiateCommand(Map(task.declarations.head -> WdlString("world")), NoFunctions)
      command shouldBe Success("echo hello world")
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
      
      val ns = WdlNamespace.loadUsingSource(wdl, None, None)
      val task = ns.findTask("t").get
      val command = task.instantiateCommand(Map(task.declarations.head -> WdlOptionalValue.none(WdlStringType)), NoFunctions)
      command.get shouldBe "echo"
    }
  }
}
