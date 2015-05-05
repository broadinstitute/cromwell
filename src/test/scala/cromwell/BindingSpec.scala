package cromwell

import java.io.File

import cromwell.binding.WdlValue
import cromwell.binding.types.{WdlIntegerType, WdlStringType, WdlFileType}
import org.scalatest.{FlatSpec, Matchers}

class BindingSpec extends FlatSpec with Matchers with WdlThreeStepFixture {
    "Binding Workflow" should "Have correct name for workflow" in {
        binding.workflow.name shouldEqual "three_step"
    }
    it should "Have correct FQN" in {
        binding.workflow.fullyQualifiedName shouldEqual "three_step"
    }
    it should "Have no parent" in {
        binding.workflow.parent shouldEqual None
    }
    it should "Have three 'Call' objects" in {
        binding.workflow.calls.size shouldEqual 3
    }

    "Binding Tasks" should "Have three task definitions" in {
        binding.tasks.size shouldEqual 3
    }
    it should "Have a task with name 'wc'" in {
        val task = binding.findTask("wc") getOrElse fail("No 'wc' task found")
        task.name shouldEqual "wc"
        task.inputs shouldEqual Map("in_file" -> WdlFileType)
        task.command.instantiate(Map("in_file" -> new WdlValue(new File("/path/to/file"), WdlFileType))) shouldEqual "wc -l /path/to/file"
        task.outputs.size shouldEqual 1
        task.outputs.head.name shouldEqual "count"
        task.outputs.head.wdlType shouldEqual WdlIntegerType
    }
    it should "Have a task with name 'cgrep'" in {
        val task = binding.findTask("cgrep") getOrElse fail("No 'cgrep' task found")
        task.name shouldEqual "cgrep"
        task.inputs shouldEqual Map("pattern" -> WdlStringType, "in_file" -> WdlFileType)
        task.command.instantiate(
            Map("pattern" -> new WdlValue("^...$", WdlStringType), "in_file" -> new WdlValue(new File("/path/to/file"), WdlFileType))
        ) shouldEqual "grep '^...$' /path/to/file | wc -l"
        task.outputs.size shouldEqual 1
        task.outputs.head.name shouldEqual "count"
        task.outputs.head.wdlType shouldEqual WdlIntegerType
    }
    it should "Have a task with name 'ps'" in {
        val task = binding.findTask("ps") getOrElse fail("No 'ps' task found")
        task.name shouldEqual "ps"
        task.inputs shouldEqual Map()
        task.command.instantiate(Map()) shouldEqual "ps"
        task.outputs.size shouldEqual 1
        task.outputs.head.name shouldEqual "procs"
        task.outputs.head.wdlType shouldEqual WdlFileType
    }
}
