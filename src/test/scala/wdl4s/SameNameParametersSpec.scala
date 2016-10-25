package wdl4s

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.expression.NoFunctions
import wdl4s.values.WdlString

class SameNameParametersSpec extends FlatSpec with Matchers {
  val namespace1 = WdlNamespaceWithWorkflow.load(
    """
       |task test {
       |  String x
       |  command { ./script ${x} ${x} ${x} }
       |}
       |workflow wf { call test }
     """.stripMargin
  )
  val task = namespace1.findTask("test").get

  "A task with command that uses the same parameter more than once" should "only count it as one input" in {
    task.declarations.map(_.toWdlString) shouldEqual Seq("String x")
  }

  it should "instantiate the command with duplicated parameter names properly" in {
    val inputs = task.inputsFromMap(Map("test.x" -> WdlString("foo")))
    task.instantiateCommand(inputs, NoFunctions).get shouldEqual "./script foo foo foo"
  }
}
