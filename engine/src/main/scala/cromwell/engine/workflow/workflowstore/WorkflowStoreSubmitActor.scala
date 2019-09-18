package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import common.validation.IOChecked._
import common.validation.Validation._
import cromwell.core.Dispatcher._
import cromwell.core._
import cromwell.engine.instrumentation.WorkflowInstrumentation
import cromwell.engine.workflow.WorkflowMetadataHelper
import cromwell.engine.workflow.WorkflowProcessingEventPublishing._
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowStoreState.WorkflowStoreState
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.{WorkflowStoreState, WorkflowSubmissionResponse}
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.{WorkflowSubmitFailed, WorkflowSubmittedToStore, WorkflowsBatchSubmittedToStore}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import spray.json._

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

final case class WorkflowStoreSubmitActor(store: WorkflowStore, serviceRegistryActor: ActorRef) extends Actor 
  with ActorLogging with WorkflowMetadataHelper with MonitoringCompanionHelper with WorkflowInstrumentation {
  implicit val ec: ExecutionContext = context.dispatcher

  val workflowStoreReceive: Receive = {
    case cmd: SubmitWorkflow =>
      addWork()
      val sndr = sender()

      val futureResponse = for {
        submissionResponses <- storeWorkflowSources(NonEmptyList.of(cmd.source))
        id = submissionResponses.head.id
        _ = registerSubmission(id, cmd.source)
      } yield WorkflowSubmissionResponse(submissionResponses.head.state, id)

      futureResponse onComplete {
        case Success(workflowSubmissionResponse) =>
          val wfType = cmd.source.workflowType.getOrElse("Unspecified type")
          val wfTypeVersion = cmd.source.workflowTypeVersion.getOrElse("Unspecified version")
          log.info("{} ({}) workflow {} submitted", wfType, wfTypeVersion, workflowSubmissionResponse.id)
          val labelsMap = convertJsonToLabelsMap(cmd.source.labelsJson)
          publishLabelsToMetadata(workflowSubmissionResponse.id, labelsMap, serviceRegistryActor).toErrorOr.toTry.get
          sndr ! WorkflowSubmittedToStore(
            workflowSubmissionResponse.id,
            convertDatabaseStateToApiState(workflowSubmissionResponse.state)
          )
          removeWork()
        case Failure(throwable) =>
          log.error("Workflow {} submit failed.", throwable)
          sndr ! WorkflowSubmitFailed(throwable)
          removeWork()
      }

    case cmd: BatchSubmitWorkflows =>
      addWork()
      val sndr = sender()

      val futureResponses = for {
        submissionResponses <- storeWorkflowSources(cmd.sources)
        _ = (submissionResponses.toList.map(res => res.id) zip cmd.sources.toList) foreach (registerSubmission _).tupled
      } yield submissionResponses

      futureResponses onComplete {
        case Success(workflowSubmissionResponses) =>
          log.info("Workflows {} submitted.", workflowSubmissionResponses.toList.map(res => res.id).mkString(", "))
          val labelsMap = convertJsonToLabelsMap(cmd.sources.head.labelsJson)
          workflowSubmissionResponses.map(res => publishLabelsToMetadata(res.id, labelsMap, serviceRegistryActor))
          sndr ! WorkflowsBatchSubmittedToStore(
            workflowSubmissionResponses.map(res => res.id),
            convertDatabaseStateToApiState(workflowSubmissionResponses.head.state)
          )
          removeWork()
        case Failure(throwable) =>
          log.error("Workflow {} submit failed.", throwable)
          sndr ! WorkflowSubmitFailed(throwable)
          removeWork()
      }
  }
  
  override def receive = workflowStoreReceive.orElse(monitoringReceive)

  private def convertDatabaseStateToApiState(workflowStoreState: WorkflowStoreState): WorkflowState ={
    workflowStoreState match {
      case WorkflowStoreState.Submitted => WorkflowSubmitted
      case WorkflowStoreState.OnHold => WorkflowOnHold
      case WorkflowStoreState.Aborting => WorkflowAborting
      case WorkflowStoreState.Running => WorkflowRunning
    }
  }

  private def storeWorkflowSources(sources: NonEmptyList[WorkflowSourceFilesCollection]): Future[NonEmptyList[WorkflowSubmissionResponse]] = {
    for {
      workflowSubmissionResponses <- store.add(sources)
    } yield workflowSubmissionResponses
  }

  private def convertJsonToLabelsMap(json: String): Map[String, String] = {
    json.parseJson match {
      case JsObject(inputs) => inputs.collect({
        case (key, JsString(value)) => key -> value
      })
      case _ => Map.empty
    }
  }

  /**
    * Runs processing on workflow source files before they are stored.
    *
    * @param processOptions How to process the workflow options
    * @param source         Original workflow source
    * @return Attempted updated workflow source
    */
  private def processSource(processOptions: WorkflowOptions => WorkflowOptions)
                           (source: WorkflowSourceFilesCollection): WorkflowSourceFilesCollection = {

    source.setOptions(processOptions(source.workflowOptions))
  }

  /**
    * Takes the workflow id and sends it over to the metadata service w/ default empty values for inputs/outputs
    */
  private def registerSubmission(
      id: WorkflowId,
      originalSourceFiles: WorkflowSourceFilesCollection): Unit = {
    // Increment the workflow submitted count
    incrementWorkflowState(WorkflowSubmitted)

    val actualWorkflowState = if(originalSourceFiles.workflowOnHold)
      WorkflowOnHold
    else
      WorkflowSubmitted

    val sourceFiles = processSource(_.clearEncryptedValues)(originalSourceFiles)

      val submissionEvents: List[MetadataEvent] = List(
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionTime), MetadataValue(OffsetDateTime.now)),
        MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Inputs)),
        MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Outputs)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.Status), MetadataValue(actualWorkflowState)),

        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Workflow), MetadataValue(sourceFiles.workflowSource.orNull)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_WorkflowUrl), MetadataValue(sourceFiles.workflowUrl.orNull)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Root), MetadataValue(sourceFiles.workflowRoot.orNull)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Inputs), MetadataValue(sourceFiles.inputsJson)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Options), MetadataValue(sourceFiles.workflowOptions.asPrettyJson)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Labels), MetadataValue(sourceFiles.labelsJson))
      )

      // Don't publish metadata for either workflow type or workflow type version if not defined.
      val workflowTypeAndVersionEvents: List[Option[MetadataEvent]] = List(
        sourceFiles.workflowType map { wt => MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_WorkflowType), MetadataValue(wt)) },
        sourceFiles.workflowTypeVersion map { wtv => MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_WorkflowTypeVersion), MetadataValue(wtv)) }
      )

      serviceRegistryActor ! PutMetadataAction(submissionEvents ++ workflowTypeAndVersionEvents.flatten)
      ()
  }
}

object WorkflowStoreSubmitActor {
  def props(workflowStoreDatabase: WorkflowStore, serviceRegistryActor: ActorRef) = {
    Props(WorkflowStoreSubmitActor(workflowStoreDatabase, serviceRegistryActor)).withDispatcher(ApiDispatcher)
  }

  sealed trait WorkflowStoreSubmitActorResponse
  final case class WorkflowSubmittedToStore(workflowId: WorkflowId, state: WorkflowState) extends WorkflowStoreSubmitActorResponse
  final case class WorkflowsBatchSubmittedToStore(workflowIds: NonEmptyList[WorkflowId], state: WorkflowState) extends WorkflowStoreSubmitActorResponse

  final case class WorkflowSubmitFailed(throwable: Throwable) extends WorkflowStoreSubmitActorResponse
}
