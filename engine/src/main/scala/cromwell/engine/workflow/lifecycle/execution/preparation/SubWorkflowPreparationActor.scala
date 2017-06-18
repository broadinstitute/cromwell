package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.actor.{Actor, Props}
import cromwell.backend.BackendJobBreadCrumb
import cromwell.core.Dispatcher._
import cromwell.core.WorkflowId
import cromwell.core.logging.WorkflowLogging
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.SubWorkflowKey
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActorData
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation.{CallPreparationFailed, Start, _}
import cromwell.engine.workflow.lifecycle.execution.preparation.SubWorkflowPreparationActor.SubWorkflowPreparationSucceeded
import wdl4s._
import wdl4s.values.WdlValue

class SubWorkflowPreparationActor(executionData: WorkflowExecutionActorData,
                                   callKey: SubWorkflowKey,
                                   subWorkflowId: WorkflowId) extends Actor with WorkflowLogging {
  
  private val workflowDescriptor = executionData.workflowDescriptor
  lazy val outputStore = executionData.outputStore
  lazy val expressionLanguageFunctions = executionData.expressionLanguageFunctions

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
    case Start =>
      val evaluatedInputs = resolveAndEvaluateInputs(callKey, workflowDescriptor, expressionLanguageFunctions, outputStore)
      val response = evaluatedInputs map { prepareExecutionActor }
      context.parent ! (response recover { case f => CallPreparationFailed(callKey, f) }).get
      context stop self

    case unhandled => workflowLogger.warn(self.path.name + " received an unhandled message: " + unhandled)
  }
}

object SubWorkflowPreparationActor {
  case class SubWorkflowPreparationSucceeded(workflowDescriptor: EngineWorkflowDescriptor, inputs: EvaluatedTaskInputs) extends CallPreparationActorResponse

  def props(executionData: WorkflowExecutionActorData,
            key: SubWorkflowKey,
            subWorkflowId: WorkflowId) = {
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new SubWorkflowPreparationActor(executionData, key, subWorkflowId)).withDispatcher(EngineDispatcher)
  }
}
