package cromwell.binding

import cromwell.binding.types.WdlStringType
import cromwell.binding.values.WdlString
import org.scalatest.{FlatSpec, Matchers}

class SameNameParametersSpec extends FlatSpec with Matchers {
  val namespace1 = NamespaceWithWorkflow.load(
    """
       |task test {
       |  command { ./script ${String x} ${String x} ${x} }
       |}
       |workflow wf { call test }
     """.stripMargin
  )

  "A task with command that uses the same parameter more than once" should "only count it as one input" in {
    namespace1.findTask("test").get.inputs shouldEqual Seq(TaskInput(name="x", wdlType=WdlStringType))
  }

  it should "instantiate the command with duplicated parameter names properly" in {
    namespace1.findTask("test").get.command.instantiate(Map("x" -> WdlString("foo"))).get shouldEqual "./script foo foo foo"
  }
}
