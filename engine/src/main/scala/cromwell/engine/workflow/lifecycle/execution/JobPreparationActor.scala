package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{Actor, ActorRef, Props}
import cromwell.backend._
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.logging.WorkflowLogging
import cromwell.core.{CallKey, JobKey, WorkflowId}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.CallPreparationActor._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.SubWorkflowKey
import wdl4s._
import wdl4s.exception.VariableLookupException
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values.WdlValue

import scala.util.{Failure, Success, Try}

abstract class CallPreparationActor(val workflowDescriptor: EngineWorkflowDescriptor,
                                    val outputStore: OutputStore,
                                    callKey: CallKey) extends Actor with WorkflowLogging {
  lazy val workflowIdForLogging = workflowDescriptor.id
  def expressionLanguageFunctions: WdlStandardLibraryFunctions
  def prepareExecutionActor(inputEvaluation: Map[Declaration, WdlValue]): CallPreparationActorResponse
  
  override def receive = {
    case Start =>
      val response = resolveAndEvaluateInputs() map { prepareExecutionActor }
      context.parent ! (response recover { case f => CallPreparationFailed(callKey, f) }).get
      context stop self

    case unhandled => workflowLogger.warn(self.path.name + " received an unhandled message: " + unhandled)
  }

  def resolveAndEvaluateInputs(): Try[Map[Declaration, WdlValue]] = {
    Try {
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
      )
    } recoverWith {
      case t: Throwable => Failure(new VariableLookupException(s"Couldn't resolve all inputs for ${callKey.scope.fullyQualifiedName} at index ${callKey.index}.", List(t)))
    }
  }
}

final case class JobPreparationActor(executionData: WorkflowExecutionActorData,
                                     jobKey: BackendJobDescriptorKey,
                                     factory: BackendLifecycleActorFactory,
                                     initializationData: Option[BackendInitializationData],
                                     serviceRegistryActor: ActorRef,
                                     backendSingletonActor: Option[ActorRef])
  extends CallPreparationActor(executionData.workflowDescriptor, executionData.outputStore, jobKey) {

  override lazy val expressionLanguageFunctions = factory.expressionLanguageFunctions(workflowDescriptor.backendDescriptor, jobKey, initializationData)
  
  override def prepareExecutionActor(inputEvaluation: Map[Declaration, WdlValue]): CallPreparationActorResponse = {
    import RuntimeAttributeDefinition.{addDefaultsToAttributes, evaluateRuntimeAttributes}
    val curriedAddDefaultsToAttributes = addDefaultsToAttributes(factory.runtimeAttributeDefinitions(initializationData), workflowDescriptor.backendDescriptor.workflowOptions) _

    (for {
      unevaluatedRuntimeAttributes <- Try(jobKey.call.task.runtimeAttributes)
      evaluatedRuntimeAttributes <- evaluateRuntimeAttributes(unevaluatedRuntimeAttributes, expressionLanguageFunctions, inputEvaluation)
      attributesWithDefault = curriedAddDefaultsToAttributes(evaluatedRuntimeAttributes)
      jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, attributesWithDefault, inputEvaluation)
    } yield BackendJobPreparationSucceeded(jobDescriptor, factory.jobExecutionActorProps(jobDescriptor, initializationData, serviceRegistryActor, backendSingletonActor))) match {
      case Success(s) => s
      case Failure(f) => CallPreparationFailed(jobKey, f)
    }
  }
}

final case class SubWorkflowPreparationActor(executionData: WorkflowExecutionActorData,
                                             key: SubWorkflowKey,
                                             subWorkflowId: WorkflowId)
  extends CallPreparationActor(executionData.workflowDescriptor, executionData.outputStore, key) {

  override lazy val expressionLanguageFunctions = executionData.expressionLanguageFunctions

  override def prepareExecutionActor(inputEvaluation: Map[Declaration, WdlValue]): CallPreparationActorResponse = {
    val oldBackendDescriptor = workflowDescriptor.backendDescriptor
    
    val newBackendDescriptor = oldBackendDescriptor.copy(
      id = subWorkflowId,
      workflow = key.scope.calledWorkflow,
      knownValues = workflowDescriptor.knownValues ++ (inputEvaluation map { case (k, v) => k.fullyQualifiedName -> v }),
      breadCrumbs = oldBackendDescriptor.breadCrumbs :+ BackendJobBreadCrumb(workflowDescriptor.workflow, workflowDescriptor.id, key)
    )
    val engineDescriptor = workflowDescriptor.copy(backendDescriptor = newBackendDescriptor, parentWorkflow = Option(workflowDescriptor))
    SubWorkflowPreparationSucceeded(engineDescriptor, inputEvaluation)
  }
}

object CallPreparationActor {
  sealed trait CallPreparationActorCommands
  case object Start extends CallPreparationActorCommands

  sealed trait CallPreparationActorResponse
  
  case class BackendJobPreparationSucceeded(jobDescriptor: BackendJobDescriptor, bjeaProps: Props) extends CallPreparationActorResponse
  case class SubWorkflowPreparationSucceeded(workflowDescriptor: EngineWorkflowDescriptor, inputs: EvaluatedTaskInputs) extends CallPreparationActorResponse
  case class JobCallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse
  case class CallPreparationFailed(jobKey: JobKey, throwable: Throwable) extends CallPreparationActorResponse
}

object JobPreparationActor {
  def props(executionData: WorkflowExecutionActorData,
            jobKey: BackendJobDescriptorKey,
            factory: BackendLifecycleActorFactory,
            initializationData: Option[BackendInitializationData],
            serviceRegistryActor: ActorRef,
            backendSingletonActor: Option[ActorRef]) = {
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new JobPreparationActor(executionData,
      jobKey,
      factory,
      initializationData,
      serviceRegistryActor,
      backendSingletonActor)).withDispatcher(EngineDispatcher)
  }
}

object SubWorkflowPreparationActor {
  def props(executionData: WorkflowExecutionActorData,
            key: SubWorkflowKey,
            subWorkflowId: WorkflowId) = {
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new SubWorkflowPreparationActor(executionData, key, subWorkflowId)).withDispatcher(EngineDispatcher)
  }
}
