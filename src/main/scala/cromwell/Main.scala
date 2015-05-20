package cromwell

import cromwell.server.CromwellServer

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.concurrent.ExecutionContext.Implicits.global

import java.io.File
import java.nio.file.Paths
import akka.pattern.ask
import akka.actor.ActorSystem
import akka.util.Timeout
import cromwell.binding._
import cromwell.binding.values.{WdlString, WdlFile, WdlValue}
import cromwell.engine.WorkflowActor.GetOutputs
import cromwell.engine.{ActorWorkflowManager, WorkflowId, WorkflowManagerActor}
import cromwell.engine.WorkflowManagerActor.{WorkflowStatus, WorkflowOutputs, SubmitWorkflow}

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
    getAction(args.headOption) match {
      case Some(x) if x == Actions.parse => parseIt(args.tail)
      case Some(x) if x == Actions.run => runIt(args.tail)
      case Some(x) if x == Actions.server => CromwellServer // Mention it so it gets instantiated
      case _ => usageAndExit()
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

  /** YAP */
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

  /** YAP */
//  def serverIt(args: Array[String]): Unit = {
//    val systemName = "cromwell-system"
//    val actorSystem = ActorSystem(systemName)
//
//    actorSystem.registerOnTermination({
//      actorSystem.log.info(s"$systemName shutting down")
//    })
//
//    val workflowManagerActor = actorSystem.actorOf(ActorWorkflowManager.props)
//
//    // TODO: Awaiting spray stuff, and really - once that's here this whole thing should move somewhere else
//
//
//    println("Cromwell server started")
//    // FIXME: Just testing for now
//    val HelloWdl =
//      """
//        |task hello {
//        |  command {
//        |    echo "Hello ${addressee}!"
//        |  }
//        |  output {
//        |    String salutation = read_string("stdout")
//        |  }
//        |}
//        |
//        |workflow hello {
//        |  call hello
//        |}
//      """.stripMargin
//
//    val Addressee = "hello.hello.addressee"
//    val HelloInputs: Map[FullyQualifiedName, WdlValue] = Map(Addressee -> WdlString("world"))
//    implicit val timeout = Timeout(30.seconds)
//
//    println(s"Manager is: $workflowManagerActor")
//
//    // FIXME: If you comment on the below in the PR I will end you
//    val foo = workflowManagerActor ? SubmitWorkflow(HelloWdl, HelloInputs)
//    val id = Await.result(foo.mapTo[Try[WorkflowId]], 5.seconds).get
//    println(s"ID IS $id")
//
//    val s1 = workflowManagerActor ? WorkflowStatus(id)
//    s1 foreach println
//
//    Thread.sleep(3000)
//    val z = workflowManagerActor ? WorkflowOutputs(id)
//    z.mapTo[Option[Map[String, WdlValue]]] onComplete {
//      case Success(s) =>
//        println(s)
//        val s2 = workflowManagerActor ? WorkflowStatus(id)
//        s2 foreach println
//      case Failure(e) => throw e
//    }
//  }

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
