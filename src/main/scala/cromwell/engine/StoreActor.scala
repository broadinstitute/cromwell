package cromwell.engine

import akka.actor.{Actor, Props}
import akka.event.{Logging, LoggingReceive}
import akka.pattern.pipe
import akka.util.Timeout
import cromwell.binding.values.WdlValue
import cromwell.binding.{WdlBinding, _}
import cromwell.engine.StoreActor._
import cromwell.parser.WdlParser.{Ast, Terminal}
import cromwell.util.TryUtil

import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.Try


object StoreActor {
  def props(binding: WdlBinding, inputs: WorkflowCoercedInputs) =
    Props(new StoreActor(binding, inputs))

  sealed trait StoreActorMessage
  case class CallCompleted(call: Call, callOutputs: Map[String, WdlValue]) extends StoreActorMessage
  case object FindRunnableCalls extends StoreActorMessage
  case class PrepareToStartCalls(runnableCalls: Iterable[Call]) extends StoreActorMessage
  case class UpdateStatus(call: Call, status: ExecutionStatus.Value) extends StoreActorMessage
  case object GetOutputs extends StoreActorMessage
  case class GetLocallyQualifiedInputs(call: Call) extends StoreActorMessage
  implicit val ActorTimeout = Timeout(5 seconds)
}

/** Actor to hold symbol and execution status data for a single workflow. */
class StoreActor(binding: WdlBinding, inputs: WorkflowCoercedInputs) extends Actor {
  private val symbolStore = new SymbolStore(binding, inputs)
  private val executionStore = new ExecutionStore(binding)
  private val log = Logging(context.system, this)

  override def receive: Receive = LoggingReceive {

    case CallCompleted(call, callOutputs) =>
      sender ! handleCallCompleted(call, callOutputs)

    case FindRunnableCalls =>
      sender ! executionStore.runnableCalls

    case PrepareToStartCalls(calls) =>
      handlePrepareToStartCalls(calls) pipeTo sender

    case GetOutputs =>
      sender ! (symbolStore.getOutputs map symbolStoreEntryToMapEntry).toMap

    case GetLocallyQualifiedInputs(call) =>
      sender ! symbolStore.locallyQualifiedInputs(call)

    case UpdateStatus(call, status) =>
      executionStore.updateStatus(call, status)
  }

  private def symbolStoreEntryToMapEntry(e: SymbolStoreEntry): (String, WdlValue) =
    e.key.scope + "." + e.key.name -> e.wdlValue.get

  private def updateOutputs(call: Call, callOutputs: Map[String, WdlValue]): Unit = {

    def addOutputValueToSymbolStore(callOutput: (String, WdlValue)): Try[Unit] =
      symbolStore.addOutputValue(call.fullyQualifiedName, callOutput._1, Some(callOutput._2), callOutput._2.wdlType)

    val addedEntries = callOutputs map addOutputValueToSymbolStore
    val failureMessages = TryUtil.stringifyFailures(addedEntries)

    if (failureMessages.isEmpty) {
      executionStore.updateStatus(call, ExecutionStatus.Done)
    } else {
      val errorMessages = failureMessages.mkString("\n")
      log.error(errorMessages)
      executionStore.updateStatus(call, ExecutionStatus.Failed)
      throw new RuntimeException(errorMessages)
    }
  }

  /**
   * Updates outputs for the completed call and returns a `Future[Boolean]` which is true
   * if the workflow is now "done".  Current "done" for a workflow means all calls are done.
   */
  private def handleCallCompleted(call: Call, callOutputs: Map[String, WdlValue]): Boolean = {
    updateOutputs(call, callOutputs)
    executionStore.isWorkflowDone
  }

  private def handlePrepareToStartCalls(calls: Iterable[Call]): Future[Unit] = Future {
    def copyOutputsToInputs(call: Call) = {
      def copyOutputToInput(inputName: String, expression: WdlExpression) = {
        val ast = expression.ast.asInstanceOf[Ast]
        val Seq(lhs, rhs) = Seq("lhs", "rhs").map {
          ast.getAttribute(_).asInstanceOf[Terminal].getSourceString
        }
        val outputFqn = Seq(call.parent.get.name, lhs, rhs).mkString(".")
        val inputFqn = Seq(call.parent.get.name, call.name, inputName).mkString(".")
        symbolStore.copyOutputToInput(outputFqn, inputFqn)
      }

      call.inputMappings.foreach { case (inputName, expression) =>
        copyOutputToInput(inputName, expression)
      }
    }
    calls.foreach { call => copyOutputsToInputs(call) }
  }
}
