package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.Props
import cromwell.backend.BackendJobDescriptor
import cromwell.core.{CallKey, JobKey}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.OutputStore
import wdl4s.wdl.exception.VariableLookupException
import wdl4s.wdl.expression.WdlStandardLibraryFunctions
import wdl4s.wdl.values.WdlValue
import wdl4s.wdl.{Declaration, Scatter}

import scala.util.{Failure, Try}

object CallPreparation {
  sealed trait CallPreparationActorCommands
  case object Start extends CallPreparationActorCommands

  trait CallPreparationActorResponse

  case class BackendJobPreparationSucceeded(jobDescriptor: BackendJobDescriptor, bjeaProps: Props) extends CallPreparationActorResponse

  case class JobCallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse
  case class CallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse

  def resolveAndEvaluateInputs(callKey: CallKey,
                               workflowDescriptor: EngineWorkflowDescriptor,
                               expressionLanguageFunctions: WdlStandardLibraryFunctions,
                               outputStore: OutputStore): Try[Map[Declaration, WdlValue]] = {
    val call = callKey.scope
    val scatterMap = callKey.index flatMap { i =>
      // Will need update for nested scatters
      call.ancestry collectFirst { case s: Scatter => Map(s -> i) }
    } getOrElse Map.empty[Scatter, Int]

    call.evaluateTaskInputs(
      workflowDescriptor.backendDescriptor.knownValues,
      expressionLanguageFunctions,
      outputStore.fetchNodeOutputEntries,
      scatterMap
    ) recoverWith {
      case t: Throwable => Failure[Map[Declaration, WdlValue]](new VariableLookupException(s"Couldn't resolve all inputs for ${callKey.scope.fullyQualifiedName} at index ${callKey.index}.", List(t)))
    }
  }
}
