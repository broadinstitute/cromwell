package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.{Actor, Props}
import cats.data.Validated.{Invalid, Valid}
import cromwell.backend.BackendJobBreadCrumb
import cromwell.core.Dispatcher._
import cromwell.core.WorkflowId
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.SubWorkflowKey
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation.{CallPreparationFailed, Start, _}
import cromwell.engine.workflow.lifecycle.execution.preparation.SubWorkflowPreparationActor.SubWorkflowPreparationSucceeded
import cromwell.engine.{EngineWorkflowDescriptor, WdlFunctions}
import lenthall.exception.MessageAggregation
import wom.values.WomEvaluatedCallInputs

class SubWorkflowPreparationActor(workflowDescriptor: EngineWorkflowDescriptor,
                                  expressionLanguageFunctions: WdlFunctions,
                                  callKey: SubWorkflowKey,
                                  subWorkflowId: WorkflowId) extends Actor with WorkflowLogging {

  lazy val workflowIdForLogging = workflowDescriptor.id

  def prepareExecutionActor(inputEvaluation: WomEvaluatedCallInputs): CallPreparationActorResponse = {
    val oldBackendDescriptor = workflowDescriptor.backendDescriptor

    val newBackendDescriptor = oldBackendDescriptor.copy(
      id = subWorkflowId,
      workflow = callKey.node.callable,
      // TODO WOM: need FQN for input definitions somehow ? For now don't worry about subWF
      knownValues = Map.empty,//workflowDescriptor.knownValues ++ (inputEvaluation map { case (k, v) => k.fullyQualifiedName -> v }),
      breadCrumbs = oldBackendDescriptor.breadCrumbs :+ BackendJobBreadCrumb(workflowDescriptor.workflow, workflowDescriptor.id, callKey)
    )
    val engineDescriptor = workflowDescriptor.copy(backendDescriptor = newBackendDescriptor, parentWorkflow = Option(workflowDescriptor))
    SubWorkflowPreparationSucceeded(engineDescriptor, inputEvaluation)
  }

  override def receive = {
    case Start(outputStore) =>
      val evaluatedInputs = resolveAndEvaluateInputs(callKey, workflowDescriptor, expressionLanguageFunctions, outputStore)
      evaluatedInputs map { prepareExecutionActor } match {
        case Valid(response) => context.parent ! response
        case Invalid(f) => context.parent ! CallPreparationFailed(callKey, new MessageAggregation {
          override def exceptionContext: String = "Failed to evaluate inputs for sub workflow"
          override def errorMessages: Traversable[String] = f.toList
        })
      }
      context stop self

    case unhandled => workflowLogger.warn(self.path.name + " received an unhandled message: " + unhandled)
  }
}

object SubWorkflowPreparationActor {
  case class SubWorkflowPreparationSucceeded(workflowDescriptor: EngineWorkflowDescriptor, inputs: WomEvaluatedCallInputs) extends CallPreparationActorResponse

  def props(workflowDescriptor: EngineWorkflowDescriptor,
            expressionLanguageFunctions: WdlFunctions,
            key: SubWorkflowKey,
            subWorkflowId: WorkflowId) = {
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new SubWorkflowPreparationActor(workflowDescriptor, expressionLanguageFunctions, key, subWorkflowId)).withDispatcher(EngineDispatcher)
  }
}
