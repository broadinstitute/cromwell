package wdl

import common.validation.Validation._
import org.scalatest.{FlatSpec, Matchers}
import wdl.expression.NoFunctions
import wom.values.WomString

class AdvancedInterpolationSpec extends FlatSpec with Matchers {
  val wdl =
    """
      |task test {
      |  String eval_this
      |  String var="inside"
      |  String evaled="${eval_this}"
      |
      |  command {
      |    echo '${eval_this} ${evaled} ${var}'
      |  }
      |  output {
      |    String out = read_string(stdout())
      |  }
      |}
      |
      |workflow testWF {
      |  String var = "outside"
      |  call test
      |  call test as test2 {input: eval_this="${var}"}
      |}
    """.stripMargin

  val namespace = WdlNamespaceWithWorkflow.load(wdl, Seq.empty).get

  it should "be able to resolve interpolated strings within interpolated strings" in {
    val testCall = namespace.workflow.taskCalls.find(_.unqualifiedName == "test") getOrElse {
      fail("call 'test' not found")
    }
    val inputs = testCall.evaluateTaskInputs(Map("testWF.test.eval_this" -> WomString("${var}")), NoFunctions)
    val testCmd = testCall.task.instantiateCommand(inputs.get, NoFunctions).toTry.get.head
    testCmd.commandString shouldEqual "echo 'inside inside inside'"
  }
}
