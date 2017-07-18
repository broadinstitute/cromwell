package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._
import cromwell.core.Dispatcher._
import cromwell.core._
import cromwell.engine.workflow.lifecycle.execution.WorkflowMetadataHelper
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.{WorkflowSubmitFailed, WorkflowSubmittedToStore, WorkflowsBatchSubmittedToStore}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

final case class WorkflowStoreSubmitActor(store: WorkflowStore, serviceRegistryActor: ActorRef) extends Actor with ActorLogging with WorkflowMetadataHelper with MonitoringCompanionHelper {
  implicit val ec: ExecutionContext = context.dispatcher

  val workflowStoreReceive: Receive = {
    case cmd: SubmitWorkflow =>
      addWork()
      val sndr = sender()

      val futureId = for {
        ids <- storeWorkflowSources(NonEmptyList.of(cmd.source))
        id = ids.head
        _ <- registerSubmissionWithMetadataService(id, cmd.source)
      } yield id

      futureId onComplete {
        case Success(id) =>
          log.info("Workflow {} submitted.", id)
          sndr ! WorkflowSubmittedToStore(id)
          removeWork()
        case Failure(throwable) =>
          log.error("Workflow {} submit failed.", throwable)
          sndr ! WorkflowSubmitFailed(throwable)
          removeWork()
      }

    case cmd: BatchSubmitWorkflows =>
      addWork()
      val sndr = sender()

      val futureIds = for {
        ids <- storeWorkflowSources(cmd.sources)
        _ <- (ids.toList zip cmd.sources.toList) traverse (registerSubmissionWithMetadataService _).tupled
      } yield ids

      futureIds onComplete {
        case Success(ids) =>
          log.info("Workflows {} submitted.", ids.toList.mkString(", "))
          sndr ! WorkflowsBatchSubmittedToStore(ids)
          removeWork()
        case Failure(throwable) =>
          log.error("Workflow {} submit failed.", throwable)
          sndr ! WorkflowSubmitFailed(throwable)
          removeWork()
      }
  }
  
  override def receive = workflowStoreReceive.orElse(monitoringReceive)

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
  private def registerSubmissionWithMetadataService(
      id: WorkflowId,
      originalSourceFiles: WorkflowSourceFilesCollection): Future[Unit] = {
    processSource(_.clearEncryptedValues)(originalSourceFiles) map { sourceFiles =>
      val submissionEvents: List[MetadataEvent] = List(
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionTime), MetadataValue(OffsetDateTime.now.toString)),
        MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Inputs)),
        MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Outputs)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.Status), MetadataValue(WorkflowSubmitted)),

        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Workflow), MetadataValue(sourceFiles.workflowSource)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Inputs), MetadataValue(sourceFiles.inputsJson)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Options), MetadataValue(sourceFiles.workflowOptionsJson)),
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
}

object WorkflowStoreSubmitActor {
  def props(workflowStoreDatabase: WorkflowStore, serviceRegistryActor: ActorRef) = {
    Props(WorkflowStoreSubmitActor(workflowStoreDatabase, serviceRegistryActor)).withDispatcher(ApiDispatcher)
  }

  sealed trait WorkflowStoreSubmitActorResponse
  final case class WorkflowSubmittedToStore(workflowId: WorkflowId) extends WorkflowStoreSubmitActorResponse
  final case class WorkflowsBatchSubmittedToStore(workflowIds: NonEmptyList[WorkflowId]) extends WorkflowStoreSubmitActorResponse

  final case class WorkflowSubmitFailed(throwable: Throwable) extends WorkflowStoreSubmitActorResponse
}
