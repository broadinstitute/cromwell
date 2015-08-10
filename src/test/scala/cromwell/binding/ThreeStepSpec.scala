package cromwell.binding

import cromwell.binding.types.{WdlFileType, WdlIntegerType, WdlStringType}
import cromwell.binding.values.{WdlFile, WdlString}
import cromwell.parser.BackendType
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpec, Matchers}

class ThreeStepSpec extends FlatSpec with Matchers {
  val namespace = NamespaceWithWorkflow.load(SampleWdl.ThreeStep.wdlSource(), BackendType.LOCAL)

  "Binding Workflow" should "Have correct name for workflow" in {
    namespace.workflow.name shouldEqual "three_step"
  }
  it should "Have correct FQN" in {
    namespace.workflow.fullyQualifiedName shouldEqual "three_step"
  }
  it should "Have no parent" in {
    namespace.workflow.parent shouldEqual None
  }
  it should "Have three 'Call' objects" in {
    namespace.workflow.calls.size shouldEqual 3
  }
  it should "Have three outputs" in {
    namespace.workflow.outputs.size shouldEqual 3
  }

  "Binding Tasks" should "Have three task definitions" in {
    namespace.tasks.size shouldEqual 3
  }

  it should "Have a task with name 'wc'" in {
    val task = namespace.findTask("wc") getOrElse fail("No 'wc' task found")
    task.name shouldEqual "wc"
    task.inputs shouldEqual Vector(TaskInput("in_file", WdlFileType, postfixQuantifier=None))
    task.instantiateCommand(Map("in_file" -> WdlFile("/path/to/file"))).get shouldEqual "cat /path/to/file | wc -l"
    task.outputs.size shouldEqual 1
    task.outputs.head.name shouldEqual "count"
    task.outputs.head.wdlType shouldEqual WdlIntegerType
  }
  it should "Have a task with name 'cgrep'" in {
    val task = namespace.findTask("cgrep") getOrElse fail("No 'cgrep' task found")
    task.name shouldEqual "cgrep"
    task.inputs shouldEqual Vector(
      TaskInput("pattern", WdlStringType, postfixQuantifier=None),
      TaskInput("in_file", WdlFileType, postfixQuantifier=None)
    )
    task.instantiateCommand(
      Map("pattern" -> WdlString("^...$"), "in_file" -> WdlFile("/path/to/file"))
    ).get shouldEqual "grep '^...$' /path/to/file | wc -l"
    task.outputs.size shouldEqual 1
    task.outputs.head.name shouldEqual "count"
    task.outputs.head.wdlType shouldEqual WdlIntegerType
  }
  it should "Have a task with name 'ps'" in {
    val task = namespace.findTask("ps") getOrElse fail("No 'ps' task found")
    task.name shouldEqual "ps"
    task.inputs shouldEqual Vector()
    task.instantiateCommand(Map()).get shouldEqual "ps"
    task.outputs.size shouldEqual 1
    task.outputs.head.name shouldEqual "procs"
    task.outputs.head.wdlType shouldEqual WdlFileType
  }
}
