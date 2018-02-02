package wdl

import common.validation.Validation._
import org.scalatest.{FlatSpec, Matchers}
import wdl.expression.NoFunctions
import wom.values.WomString

class SameNameParametersSpec extends FlatSpec with Matchers {
  val namespace1 = WdlNamespaceWithWorkflow.load(
    """
       |task test {
       |  String x
       |  command { ./script ${x} ${x} ${x} }
       |}
       |workflow wf { call test }
     """.stripMargin, Seq.empty
  ).get
  val task = namespace1.findTask("test").get

  "A task with command that uses the same parameter more than once" should "only count it as one input" in {
    task.declarations.map(_.toWdlString) shouldEqual Seq("String x")
  }

  it should "instantiate the command with duplicated parameter names properly" in {
    val inputs = task.inputsFromMap(Map("test.x" -> WomString("foo")))
    task.instantiateCommand(inputs, NoFunctions).toTry.get.head.commandString shouldEqual "./script foo foo foo"
  }
}
