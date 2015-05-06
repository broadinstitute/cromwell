package cromwell

import java.io.File

import cromwell.binding.types.{WdlFileType, WdlStringType}
import cromwell.binding.{WdlBinding, WdlValue}
import cromwell.parser.WdlParser.SyntaxError

object Main {
  def main(args: Array[String]) {
    val actions = List("parse", "run")
    if (args.length < 2 || !actions.contains(args(0))) {
      println("Usage: cromwell.jar parse <wdl file>")
      println("Usage: cromwell.jar run <wdl file>")
      System.exit(-1)
    }
    if (args(0).equals("parse")) {
      val ast = WdlBinding.getAst(new File(args(1)))
      println(ast.toPrettyString)
    }
    if (args(0).equals("run")) {
      try {
        val binding = WdlBinding.process(new File(args(1)))
        println(s"Workflow FQN: ${binding.workflow.fullyQualifiedName}")
        println(s"Workflow Inputs: ${binding.workflow.inputs}")
        println(s"Workflow Outputs: ${binding.workflow.outputs}")
        println("")

        val params = Map(
          "in_file" -> new WdlValue(new File("/usr/share/dict/words"), WdlFileType),
          "pattern" -> new WdlValue("^...$", WdlStringType)
        )

        binding.workflow.calls.foreach { call =>
          println(s"Call: $call")
          println(s"FQN: ${call.fullyQualifiedName}")
          println(s"Unsatisfied Inputs: ${call.unsatisfiedInputs}")
          println(s"Parent: ${call.parent}")
          println("")
        }

        binding.tasks.foreach { task =>
          println(s"Task: ${task.name}")
          println(s"Abstract Command: ${task.command}")
          println(s"Inputs: ${task.command.inputs}")
          println(s"Instantiated: ${task.command.instantiate(params)}")
          println("")
        }
      } catch {
        case e: SyntaxError =>
          println(e);
          System.exit(-1);
      }
    }
  }
}
