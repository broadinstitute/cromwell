package cromwell.binding

import cromwell.binding.types.{WdlFileType, WdlIntegerType, WdlStringType}
import cromwell.binding.values.{WdlFile, WdlString}
import org.scalatest.{FlatSpec, Matchers}

class ThreeStepSpec extends FlatSpec with Matchers with ThreeStepFixture {
  "Binding Workflow" should "Have one workflow definition" in {
    binding.workflows.size shouldEqual 1
  }
  it should "Have zero imported workflow definition" in {
    binding.importedWorkflows.size shouldEqual 0
  }
  it should "Have correct name for workflow" in {
    binding.workflows.head.name shouldEqual "three_step"
  }
  it should "Have correct FQN" in {
    binding.workflows.head.fullyQualifiedName shouldEqual "three_step"
  }
  it should "Have no parent" in {
    binding.workflows.head.parent shouldEqual None
  }
  it should "Have three 'Call' objects" in {
    binding.workflows.head.calls.size shouldEqual 3
  }
  it should "Have three outputs" in {
    binding.workflows.head.outputs.size shouldEqual 3
  }

  "Binding Tasks" should "Have three task definitions" in {
    binding.tasks.size shouldEqual 3
  }
  it should "Have zero imported tasks" in {
    binding.importedTasks.size shouldEqual 0
  }
  it should "Have a task with name 'wc'" in {
    val task = binding.findTask("wc") getOrElse fail("No 'wc' task found")
    task.name shouldEqual "wc"
    task.inputs shouldEqual Map("in_file" -> WdlFileType)
    task.command.instantiate(Map("in_file" -> WdlFile("/path/to/file"))).get shouldEqual "cat /path/to/file | wc -l"
    task.outputs.size shouldEqual 1
    task.outputs.head.name shouldEqual "count"
    task.outputs.head.wdlType shouldEqual WdlIntegerType
  }
  it should "Have a task with name 'cgrep'" in {
    val task = binding.findTask("cgrep") getOrElse fail("No 'cgrep' task found")
    task.name shouldEqual "cgrep"
    task.inputs shouldEqual Map("pattern" -> WdlStringType, "in_file" -> WdlFileType)
    task.command.instantiate(
      Map("pattern" -> WdlString("^...$"), "in_file" -> WdlFile("/path/to/file"))
    ).get shouldEqual "grep '^...$' /path/to/file | wc -l"
    task.outputs.size shouldEqual 1
    task.outputs.head.name shouldEqual "count"
    task.outputs.head.wdlType shouldEqual WdlIntegerType
  }
  it should "Have a task with name 'ps'" in {
    val task = binding.findTask("ps") getOrElse fail("No 'ps' task found")
    task.name shouldEqual "ps"
    task.inputs shouldEqual Map()
    task.command.instantiate(Map()).get shouldEqual "ps"
    task.outputs.size shouldEqual 1
    task.outputs.head.name shouldEqual "procs"
    task.outputs.head.wdlType shouldEqual WdlFileType
  }
}
