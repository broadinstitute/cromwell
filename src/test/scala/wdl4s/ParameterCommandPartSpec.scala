package wdl4s

import better.files._
import wdl4s.command.ParameterCommandPart
import wdl4s.parser.WdlParser.SyntaxError

import scala.util.{Failure, Try}

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
        WdlNamespace.load(
          """task param_test {
            |  Boolean f
            |
            |  command <<<
            |  ./binary ${true="--true" f}
            |  >>>
            |}
            |
            |workflow wf {call param_test}
          """.stripMargin)
      ) match {
        case Failure(s: SyntaxError) => // expected
        case _ => fail("Expecting a syntax error")
      }
    }
  }
}
