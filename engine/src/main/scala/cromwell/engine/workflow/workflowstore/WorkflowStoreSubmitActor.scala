package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.Monad
import cats.data.EitherT._
import cats.data.NonEmptyList
import cats.instances.future._
import cats.instances.list._
import cats.instances.vector._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Parse._
import cromwell.core.Dispatcher._
import cromwell.core._
import cromwell.core.labels.{Label, Labels}
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState
import cromwell.database.sql.tables.WorkflowStoreEntry.WorkflowStoreState.WorkflowStoreState
import cromwell.engine.instrumentation.WorkflowInstrumentation
import cromwell.engine.workflow.WorkflowMetadataHelper
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowSubmissionResponse
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor._
import cromwell.engine.workflow.workflowstore.WorkflowStoreSubmitActor.{WorkflowSubmitFailed, WorkflowSubmittedToStore, WorkflowsBatchSubmittedToStore}
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import spray.json.{JsObject, JsString, JsValue, _}
import wom.core.WorkflowOptionsJson

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

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
        _ <- registerSubmission(id, cmd.source)
      } yield WorkflowSubmissionResponse(submissionResponses.head.state, id)

      futureResponse onComplete {
        case Success(futureResponse) =>
          val wfType = cmd.source.workflowType.getOrElse("Unspecified type")
          val wfTypeVersion = cmd.source.workflowTypeVersion.getOrElse("Unspecified version")
          log.info("{} ({}) workflow {} submitted", wfType, wfTypeVersion, futureResponse.id)
          //********** trying to test something **** remove before creating pull request ***********
          validateLabels(cmd.source.labelsJson) foreach { labels => publishLabelsToMetadata(futureResponse.id, labels); () }
          sndr ! WorkflowSubmittedToStore(futureResponse.id, convertDatabaseStateToApiState(futureResponse.state))
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
        _ <- (submissionResponses.toList.map(res => res.id) zip cmd.sources.toList) traverse (registerSubmission _).tupled
      } yield submissionResponses

      futureResponses onComplete {
        case Success(futureResponses) =>
          log.info("Workflows {} submitted.", futureResponses.toList.map(res => res.id).mkString(", "))
          sndr ! WorkflowsBatchSubmittedToStore(futureResponses.map(res => res.id), convertDatabaseStateToApiState(futureResponses.head.state))
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
      processedSources <- processSources(sources, _.asPrettyJson)
      workflowSubmissionResponses <- store.add(processedSources)
    } yield workflowSubmissionResponses
  }

  private def processSources(sources: NonEmptyList[WorkflowSourceFilesCollection],
                             processOptions: WorkflowOptions => WorkflowOptionsJson): Future[NonEmptyList[WorkflowSourceFilesCollection]] = {
    val nelFutures: NonEmptyList[Future[WorkflowSourceFilesCollection]] = sources map processSource(processOptions)
    val listFutures: List[Future[WorkflowSourceFilesCollection]] = nelFutures.toList
    val futureLists: Future[List[WorkflowSourceFilesCollection]] = Future.sequence(listFutures)
    futureLists.map(seq => NonEmptyList.fromList(seq).get)
  }

  private def validateLabels(json: String): ErrorOr[Labels] = {

    def toLabels(inputs: Map[String, JsValue]): ErrorOr[Labels] = {
      val vectorOfValidatedLabel: Vector[ErrorOr[Label]] = inputs.toVector map {
        case (key, JsString(s)) => Label.validateLabel(key, s)
        case (key, other) => s"Invalid label $key: $other : Labels must be strings. ${Label.LabelExpectationsMessage}".invalidNel
      }

      vectorOfValidatedLabel.sequence[ErrorOr, Label] map { validatedVectorofLabel => Labels(validatedVectorofLabel) }
    }

    Try(json.parseJson) match {
      case Success(JsObject(inputs)) => toLabels(inputs)
      case Failure(reason: Throwable) => s"Workflow contains invalid labels JSON: ${reason.getMessage}".invalidNel
      case _ => """Invalid workflow labels JSON. Expected a JsObject of "labelKey": "labelValue" values.""".invalidNel
    }
  }

  private def publishLabelsToMetadata(rootWorkflowId: WorkflowId, labels: Labels): Parse[Unit] = {
    val defaultLabel = "cromwell-workflow-id" -> s"cromwell-$rootWorkflowId"
    val customLabels = labels.asMap
    Monad[Parse].pure(labelsToMetadata(customLabels + defaultLabel, rootWorkflowId))
  }

  protected def labelsToMetadata(labels: Map[String, String], workflowId: WorkflowId): Unit = {
    labels foreach { case (k, v) =>
      serviceRegistryActor ! PutMetadataAction(MetadataEvent(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Labels}:$k"), MetadataValue(v)))
    }
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
  private def registerSubmission(
      id: WorkflowId,
      originalSourceFiles: WorkflowSourceFilesCollection): Future[Unit] = {
    // Increment the workflow submitted count
    incrementWorkflowState(WorkflowSubmitted)

    val actualWorkflowState = if(originalSourceFiles.workflowOnHold)
      WorkflowOnHold
    else
      WorkflowSubmitted

    processSource(_.clearEncryptedValues)(originalSourceFiles) map { sourceFiles =>
      val submissionEvents: List[MetadataEvent] = List(
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionTime), MetadataValue(OffsetDateTime.now.toString)),
        MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Inputs)),
        MetadataEvent.empty(MetadataKey(id, None, WorkflowMetadataKeys.Outputs)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.Status), MetadataValue(actualWorkflowState)),

        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Workflow), MetadataValue(sourceFiles.workflowSource)),
        MetadataEvent(MetadataKey(id, None, WorkflowMetadataKeys.SubmissionSection, WorkflowMetadataKeys.SubmissionSection_Root), MetadataValue(sourceFiles.workflowRoot)),
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
  final case class WorkflowSubmittedToStore(workflowId: WorkflowId, state: WorkflowState) extends WorkflowStoreSubmitActorResponse
  final case class WorkflowsBatchSubmittedToStore(workflowIds: NonEmptyList[WorkflowId], state: WorkflowState) extends WorkflowStoreSubmitActorResponse

  final case class WorkflowSubmitFailed(throwable: Throwable) extends WorkflowStoreSubmitActorResponse
}
