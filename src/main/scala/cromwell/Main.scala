package cromwell

import java.io.File
import java.nio.file.Paths

import cromwell.binding._
import cromwell.binding.values.{WdlFile, WdlString}
import cromwell.server.CromwellServer

import scala.util.{Failure, Success, Try}


object Actions extends Enumeration {
  val parse, run, server = Value
}


object Main extends App {

  getAction(args.headOption) match {
    case Some(x) if x == Actions.parse => parseIt(args.tail)
    case Some(x) if x == Actions.run => runIt(args.tail)
    case Some(x) if x == Actions.server => CromwellServer // Mention it so it gets instantiated
    case _ => usageAndExit()
  }

  def printProcessedBinding(binding: WdlBinding): Unit = {
    val workflow =
      s"""
         |Workflow FQN: ${binding.workflow.fullyQualifiedName}
         |Workflow Inputs: ${binding.workflow.inputs}
         |Workflow Outputs: ${binding.workflow.outputs}
     """.stripMargin
    println(workflow.trim + "\n")

    val params = Map(
      "in_file" -> WdlFile(Paths.get("/usr/share/dict/words")),
      "pattern" -> WdlString("^...$")
    )

    binding.workflow.calls.foreach { call =>
      val callString =
        s"""
           |Call: $call
           |FQN: ${call.fullyQualifiedName}
           |Unsatisfied Inputs: ${call.unsatisfiedInputs}
           |Parent: ${call.parent}
         """.stripMargin
      println(callString + "\n")
    }

    binding.tasks.foreach { task =>
      val taskString =
        s"""
           |Task: ${task.name}
           |Abstract Command: ${task.command}
           |Inputs: ${task.command.inputs}
           |Instantiated: ${task.command.instantiate(params)}
         """.stripMargin
      println(taskString + "\n")
    }
  }

  def parseIt(args: Array[String]): Unit = {
    if (args.isEmpty) usageAndExit()
    else println(WdlBinding.getAst(new File(args(1))).toPrettyString)
  }

  def runIt(args: Array[String]): Unit = {
    if (args.isEmpty) usageAndExit()
    else {
      Try(WdlBinding.process(new File(args(1)))) match {
        case Success(b) => printProcessedBinding(b)
        case Failure(e) =>
          println(e)
          System.exit(-1)
      }
    }
  }

  def usageAndExit(): Unit = {
    println("Usage: cromwell.jar parse <wdl file>")
    println("Usage: cromwell.jar run <wdl file>")
    println("Usage: cromwell.jar server")
    System.exit(-1)
  }

  def getAction(firstArg: Option[String]): Option[Actions.Value] = for {
    arg <- firstArg
    a <- Actions.values find { _.toString == arg }
  } yield a
}
