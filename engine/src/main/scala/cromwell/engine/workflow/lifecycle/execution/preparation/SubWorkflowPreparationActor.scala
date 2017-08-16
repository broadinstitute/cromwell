package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.{Actor, Props}
import cromwell.backend.BackendJobBreadCrumb
import cromwell.core.Dispatcher._
import cromwell.core.WorkflowId
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.SubWorkflowKey
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation.{CallPreparationFailed, Start, _}
import cromwell.engine.workflow.lifecycle.execution.preparation.SubWorkflowPreparationActor.SubWorkflowPreparationSucceeded
import cromwell.engine.{EngineWorkflowDescriptor, WdlFunctions}
import wdl4s.wdl._
import wdl4s.wdl.values.WdlValue

class SubWorkflowPreparationActor(workflowDescriptor: EngineWorkflowDescriptor,
                                  expressionLanguageFunctions: WdlFunctions,
                                  callKey: SubWorkflowKey,
                                  subWorkflowId: WorkflowId) extends Actor with WorkflowLogging {

  lazy val workflowIdForLogging = workflowDescriptor.id

  def prepareExecutionActor(inputEvaluation: Map[Declaration, WdlValue]): CallPreparationActorResponse = {
    val oldBackendDescriptor = workflowDescriptor.backendDescriptor

    val newBackendDescriptor = oldBackendDescriptor.copy(
      id = subWorkflowId,
      workflow = callKey.scope.calledWorkflow,
      knownValues = workflowDescriptor.knownValues ++ (inputEvaluation map { case (k, v) => k.fullyQualifiedName -> v }),
      breadCrumbs = oldBackendDescriptor.breadCrumbs :+ BackendJobBreadCrumb(workflowDescriptor.workflow, workflowDescriptor.id, callKey)
    )
    val engineDescriptor = workflowDescriptor.copy(backendDescriptor = newBackendDescriptor, parentWorkflow = Option(workflowDescriptor))
    SubWorkflowPreparationSucceeded(engineDescriptor, inputEvaluation)
  }

  override def receive = {
    case Start(outputStore) =>
      val evaluatedInputs = resolveAndEvaluateInputs(callKey, workflowDescriptor, expressionLanguageFunctions, outputStore)
      val response = evaluatedInputs map { prepareExecutionActor }
      context.parent ! (response recover { case f => CallPreparationFailed(callKey, f) }).get
      context stop self

    case unhandled => workflowLogger.warn(self.path.name + " received an unhandled message: " + unhandled)
  }
}

object SubWorkflowPreparationActor {
  case class SubWorkflowPreparationSucceeded(workflowDescriptor: EngineWorkflowDescriptor, inputs: EvaluatedTaskInputs) extends CallPreparationActorResponse

  def props(workflowDescriptor: EngineWorkflowDescriptor,
            expressionLanguageFunctions: WdlFunctions,
            key: SubWorkflowKey,
            subWorkflowId: WorkflowId) = {
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new SubWorkflowPreparationActor(workflowDescriptor, expressionLanguageFunctions, key, subWorkflowId)).withDispatcher(EngineDispatcher)
  }
}
