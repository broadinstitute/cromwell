package cromwell.binding

import cromwell.binding.types.WdlStringType
import cromwell.binding.values.WdlString
import cromwell.parser.BackendType
import org.scalatest.{FlatSpec, Matchers}

class SameNameParametersSpec extends FlatSpec with Matchers {
  val namespace1 = NamespaceWithWorkflow.load(
    """
       |task test {
       |  String x
       |  command { ./script ${x} ${x} ${x} }
       |}
       |workflow wf { call test }
     """.stripMargin, BackendType.LOCAL
  )
  val task = namespace1.findTask("test").get

  "A task with command that uses the same parameter more than once" should "only count it as one input" in {
    task.inputs shouldEqual Seq(TaskInput(name="x", wdlType=WdlStringType))
  }

  it should "instantiate the command with duplicated parameter names properly" in {
    task.instantiateCommand(Map("x" -> WdlString("foo"))).get shouldEqual "./script foo foo foo"
  }
}
