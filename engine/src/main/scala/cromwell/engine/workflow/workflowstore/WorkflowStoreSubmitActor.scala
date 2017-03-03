package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.core.Dispatcher._
import cromwell.core._
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.{WorkflowSubmittedToStore, WorkflowsBatchSubmittedToStore}
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.services.metadata.MetadataService.PutMetadataAction

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

final case class WorkflowStoreSubmitActor(store: WorkflowStore, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {
  implicit val ec: ExecutionContext = context.dispatcher

  override def receive = {
    case cmd: SubmitWorkflow =>
      val sndr = sender()
      storeWorkflowSources(NonEmptyList.of(cmd.source)) foreach { ids =>
        val id = ids.head
        registerSubmissionWithMetadataService(id, cmd.source)
        sndr ! WorkflowSubmittedToStore(id)
        log.info("Workflow {} submitted.", id)
      }
    case cmd: BatchSubmitWorkflows =>
      val sndr = sender()
      storeWorkflowSources(cmd.sources) foreach { ids =>
        val assignedSources = ids.toList.zip(cmd.sources.toList)
        assignedSources foreach { case (id, sourceFiles) => registerSubmissionWithMetadataService(id, sourceFiles) }
        sndr ! WorkflowsBatchSubmittedToStore(ids)
        log.info("Workflows {} submitted.", ids.toList.mkString(", "))
      }
  }

  private def storeWorkflowSources(sources: NonEmptyList[WorkflowSourceFilesCollection]): Future[NonEmptyList[WorkflowId]] = {
    for {
      processedSources <- processSources(sources, _.asPrettyJson)
      workflowIds <- store.add(processedSources)
    } yield workflowIds
  }

  private def processSources(sources: NonEmptyList[WorkflowSourceFilesCollection],
                             processOptions: WorkflowOptions => WorkflowOptionsJson): Future[NonEmptyList[WorkflowSourceFilesCollection]] = {
    val nelFutures: NonEmptyList[Future[WorkflowSourceFilesCollection]] = sources map processSource(processOptions)
    val listFutures: List[Future[WorkflowSourceFilesCollection]] = nelFutures.toList
    val futureLists: Future[List[WorkflowSourceFilesCollection]] = Future.sequence(listFutures)
    futureLists.map(seq => NonEmptyList.fromList(seq).get)
  }

  /**
    * Runs processing on workflow source files before they are stored.
    *
    * @param processOptions How to process the workflow options
    * @param source         Original workflow source
    * @return Attempted updated workflow source
    */
  private def processSource(processOptions: WorkflowOptions => WorkflowOptionsJson)
                           (source: WorkflowSourceFilesCollection): Future[WorkflowSourceFilesCollection] = {
    val options = Future {
      WorkflowOptions.fromJsonString(source.workflowOptionsJson)
    }.flatMap {
      case Success(s) => Future.successful(s)
      case Failure(regrets) => Future.failed(regrets)
    }

    options map {o => source.copyOptions(processOptions(o)) }
  }

  /**
    * Takes the workflow id and sends it over to the metadata service w/ default empty values for inputs/outputs
    */
  private def registerSubmissionWithMetadataService(id: WorkflowId, originalSourceFiles: WorkflowSourceFilesCollection): Unit = {
    processSource(_.clearEncryptedValues)(originalSourceFiles) foreach { sourceFiles =>
      val submissionEvents = List(
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionTime), MetadataValue(OffsetDateTime.now.toString)),
        MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Inputs)),
        MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Outputs)),

        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Workflow), MetadataValue(sourceFiles.wdlSource)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Inputs), MetadataValue(sourceFiles.inputsJson)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Options), MetadataValue(sourceFiles.workflowOptionsJson))
      )

      serviceRegistryActor ! PutMetadataAction(submissionEvents)
    }
  }
}

object WorkflowStoreSubmitActor {
  def props(workflowStoreDatabase: WorkflowStore, serviceRegistryActor: ActorRef) = {
    Props(WorkflowStoreSubmitActor(workflowStoreDatabase, serviceRegistryActor)).withDispatcher(ApiDispatcher)
  }

  sealed trait WorkflowStoreSubmitActorResponse
  final case class WorkflowSubmittedToStore(workflowId: WorkflowId) extends WorkflowStoreSubmitActorResponse
  final case class WorkflowsBatchSubmittedToStore(workflowIds: NonEmptyList[WorkflowId]) extends WorkflowStoreSubmitActorResponse
}

