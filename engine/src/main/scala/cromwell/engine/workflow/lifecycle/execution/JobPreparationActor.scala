package cromwell.engine.workflow.lifecycle.execution

import akka.actor.{Actor, ActorRef, Props}
import cromwell.backend._
import cromwell.core.logging.WorkflowLogging
import cromwell.core.{ExecutionStore, JobKey, OutputStore}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.workflow.lifecycle.execution.JobPreparationActor._
import wdl4s._
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values.WdlValue

import scala.util.{Failure, Success, Try}

final case class JobPreparationActor(executionData: WorkflowExecutionActorData,
                                     jobKey: BackendJobDescriptorKey,
                                     factory: BackendLifecycleActorFactory,
                                     initializationData: Option[BackendInitializationData],
                                     serviceRegistryActor: ActorRef,
                                     backendSingletonActor: Option[ActorRef])
  extends Actor with WorkflowLogging {

  lazy val workflowDescriptor: EngineWorkflowDescriptor = executionData.workflowDescriptor
  lazy val workflowId = workflowDescriptor.id
  lazy val executionStore: ExecutionStore = executionData.executionStore
  lazy val outputStore: OutputStore = executionData.outputStore
  lazy val expressionLanguageFunctions = factory.expressionLanguageFunctions(
    workflowDescriptor.backendDescriptor, jobKey, initializationData)

  override def receive = {
    case Start =>
      val response = resolveAndEvaluateInputs(jobKey, expressionLanguageFunctions) map { prepareJobExecutionActor }
      context.parent ! (response recover { case f => BackendJobPreparationFailed(jobKey, f) }).get
      context stop self

    case unhandled => workflowLogger.warn(self.path.name + " received an unhandled message: " + unhandled)
  }

  def resolveAndEvaluateInputs(jobKey: BackendJobDescriptorKey,
                               wdlFunctions: WdlStandardLibraryFunctions): Try[Map[Declaration, WdlValue]] = {
    Try {
      val call = jobKey.call
      val scatterMap = jobKey.index flatMap { i =>
        // Will need update for nested scatters
        call.upstream collectFirst { case s: Scatter => Map(s -> i) }
      } getOrElse Map.empty[Scatter, Int]

      call.evaluateTaskInputs(
        workflowDescriptor.backendDescriptor.inputs,
        expressionLanguageFunctions,
        outputStore.fetchCallOutputEntries,
        scatterMap
      )
    }
  }

  private def prepareJobExecutionActor(inputEvaluation: Map[Declaration, WdlValue]): JobPreparationActorResponse = {
    import RuntimeAttributeDefinition.{addDefaultsToAttributes, evaluateRuntimeAttributes}
    val curriedAddDefaultsToAttributes = addDefaultsToAttributes(factory.runtimeAttributeDefinitions(initializationData), workflowDescriptor.backendDescriptor.workflowOptions) _

    (for {
      unevaluatedRuntimeAttributes <- Try(jobKey.call.task.runtimeAttributes)
      evaluatedRuntimeAttributes <- evaluateRuntimeAttributes(unevaluatedRuntimeAttributes, expressionLanguageFunctions, inputEvaluation)
      attributesWithDefault = curriedAddDefaultsToAttributes(evaluatedRuntimeAttributes)
      jobDescriptor = BackendJobDescriptor(workflowDescriptor.backendDescriptor, jobKey, attributesWithDefault, inputEvaluation)
    } yield BackendJobPreparationSucceeded(jobDescriptor, factory.jobExecutionActorProps(jobDescriptor, initializationData, serviceRegistryActor, backendSingletonActor))) match {
      case Success(s) => s
      case Failure(f) => BackendJobPreparationFailed(jobKey, f)
    }
  }
}

object JobPreparationActor {
  sealed trait JobPreparationActorCommands
  case object Start extends JobPreparationActorCommands

  sealed trait JobPreparationActorResponse
  case class BackendJobPreparationSucceeded(jobDescriptor: BackendJobDescriptor, bjeaProps: Props) extends JobPreparationActorResponse
  case class BackendJobPreparationFailed(jobKey: JobKey, throwable: Throwable) extends JobPreparationActorResponse

  def props(executionData: WorkflowExecutionActorData,
            jobKey: BackendJobDescriptorKey,
            factory: BackendLifecycleActorFactory,
            initializationData: Option[BackendInitializationData],
            serviceRegistryActor: ActorRef,
            backendSingletonActor: Option[ActorRef]) = {
    // Note that JobPreparationActor doesn't run on the engine dispatcher as it mostly executes backend-side code
    // (WDL expression evaluation using Backend's expressionLanguageFunctions)
    Props(new JobPreparationActor(executionData, jobKey, factory, initializationData, serviceRegistryActor, backendSingletonActor))
  }
}
