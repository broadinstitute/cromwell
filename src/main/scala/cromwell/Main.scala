package cromwell

import java.io.File
import java.nio.file.Paths
import java.util.UUID

import akka.actor.ActorSystem
import cromwell.binding._
import cromwell.binding.values.{WdlString, WdlFile, WdlValue}
import cromwell.engine.WorkflowManagerActor
import cromwell.engine.WorkflowManagerActor.{WorkflowStatus, SubmitWorkflow}
import cromwell.server.CromwellServer
import cromwell.singleuser.SingleWorkflowRunner


import scala.util.{Failure, Success, Try}

object Actions extends Enumeration {
  val parse, run, server = Value
}

/**
 * TODO: When Hellbender's command line parser is complete, switch to using that
 *
 * FIXME: All of this should be viewed as placeholder until we get a real command line handler & have real stuff to do
 *
 * FIXME: Should change to proper logging
 */
object Main extends App {
  override def main(args: Array[String]) {
    getAction(args.headOption) match {
      case Some(x) if x == Actions.parse => parseIt(args.tail)
      case Some(x) if x == Actions.run => SingleWorkflowRunner.runIt(args.tail)
      case Some(x) if x == Actions.server => CromwellServer // Mention it so it gets loaded. Note that the actor system is also lazily initialized
      case _ => usageAndExit()
    }
  }

  /**
   * Placeholder for now
   */
  def printProcessedBinding(binding: WdlBinding): Unit = {
    println(s"Workflow FQN: ${binding.workflow.fullyQualifiedName}")
    println(s"Workflow Inputs: ${binding.workflow.inputs}")
    println(s"Workflow Outputs: ${binding.workflow.outputs}")
    println("")

    val params = Map(
      "in_file" -> WdlFile(Paths.get("/usr/share/dict/words")),
      "pattern" -> WdlString("^...$")
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
  }

  /** Placeholder */
  def parseIt(args: Array[String]): Unit = {
    if (args.isEmpty) usageAndExit()
    else println(WdlBinding.getAst(new File(args(1))).toPrettyString)
  }

  /**
   * YAP (Yet Another Placeholder)
   */
  def usageAndExit(): Unit = {
    println("Usage: cromwell.jar parse <wdl file>")
    println("Usage: cromwell.jar run <wdl file>")
    println("Usage: cromwell.jar server")
    System.exit(-1)
  }

  /** YAP */
  def getAction(firstArg: Option[String]): Option[Actions.Value] = for (arg <- firstArg; a <- Actions.values find {_.toString == arg}) yield a
}
