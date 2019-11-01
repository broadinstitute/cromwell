package cromwell.services.metadata.hybridcarbonite

import akka.actor.{ActorRef, LoggingFSM, Props, Status}
import akka.pattern.pipe
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import cromwell.core.WorkflowId
import cromwell.core.io.{AsyncIo, DefaultIoCommandBuilder}
import cromwell.services.metadata.MetadataQuery
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.hybridcarbonite.CarbonitedMetadataThawingActor._
import cromwell.services.{FailedMetadataJsonResponse, SuccessfulMetadataJsonResponse}
import cromwell.util.JsonEditor
import io.circe.parser._
import io.circe.{Json, Printer}
import spray.json._

import scala.concurrent.{ExecutionContext, Future}

final class CarbonitedMetadataThawingActor(carboniterConfig: HybridCarboniteConfig, serviceRegistry: ActorRef, ioActor: ActorRef) extends LoggingFSM[ThawingState, ThawingData] {

  implicit val ec: ExecutionContext = context.dispatcher

  val asyncIo = new AsyncIo(ioActor, DefaultIoCommandBuilder)

  startWith(PendingState, PendingData)

  when(PendingState) {
    case Event(action: BuildWorkflowMetadataJsonAction, PendingData) =>
      val gcsResponse = for {
        rawResponseFromIoActor <- asyncIo.contentAsStringAsync(carboniterConfig.makePath(action.workflowId), maxBytes = Option(30 * 1024 * 1024), failOnOverflow = true)
        parsedResponse <- Future.fromTry(parse(rawResponseFromIoActor).toTry.map(GcsMetadataResponse.apply))
      } yield parsedResponse

      gcsResponse pipeTo self

      if (action.requiresLabels) { serviceRegistry ! GetRootAndSubworkflowLabels(action.workflowId) }
      val targetState = if (action.requiresLabels) NeedsGcsAndLabelsState else NeedsOnlyGcsState

      goto(targetState) using WorkingNoData(sender(), action)
  }

  when(NeedsOnlyGcsState) {
    case Event(GcsMetadataResponse(gcsResponse), WorkingNoData(requester, action)) =>
      requester.respondWith(action, gcsResponse)
      stop()
  }

  when(NeedsGcsAndLabelsState) {
    case Event(GcsMetadataResponse(gcsResponse), WorkingWithLabelsData(requester, action, labels)) =>
      requester.respondWith(action, gcsResponse, labels)
      stop()

    case Event(RootAndSubworkflowLabelsLookupResponse(_, labels), WorkingWithGcsData(requester, action, gcsResponse)) =>
      requester.respondWith(action, gcsResponse, labels)
      stop()

    case Event(GcsMetadataResponse(gcsResponse), WorkingNoData(requester, action)) =>
      stay() using WorkingWithGcsData(requester, action, gcsResponse = gcsResponse)

    case Event(RootAndSubworkflowLabelsLookupResponse(_, labels), WorkingNoData(requester, action)) =>
      stay() using WorkingWithLabelsData(requester, action, labels = labels)

    // Failures:
    case Event(RootAndSubworkflowLabelsLookupFailed(_, failure), data: WorkingData) =>
      data.requester ! FailedMetadataJsonResponse(data.action, failure)
      stop()
  }

  whenUnhandled {
    case Event(Status.Failure(failure), data: WorkingData) =>
      data.requester ! FailedMetadataJsonResponse(data.action, failure)
      stop()
    case other =>
      log.error(s"Programmer Error: Unexpected message to ${getClass.getSimpleName} ${self.path.name}: $other")
      stay()
  }
}

object CarbonitedMetadataThawingActor {

  def props(carboniterConfig: HybridCarboniteConfig, serviceRegistry: ActorRef, ioActor: ActorRef): Props =
    Props(new CarbonitedMetadataThawingActor(carboniterConfig, serviceRegistry, ioActor))

  final case class GcsMetadataResponse(response: Json)

  sealed trait ThawingState
  case object PendingState extends ThawingState
  case object NeedsGcsAndLabelsState extends ThawingState
  case object NeedsOnlyGcsState extends ThawingState

  sealed trait ThawingData
  case object PendingData extends ThawingData
  sealed trait WorkingData extends ThawingData {
    def requester: ActorRef
    def action: BuildWorkflowMetadataJsonAction
  }
  final case class WorkingNoData(requester: ActorRef, action: BuildWorkflowMetadataJsonAction) extends WorkingData
  final case class WorkingWithLabelsData(requester: ActorRef, action: BuildWorkflowMetadataJsonAction, labels: Map[WorkflowId, Map[String, String]]) extends WorkingData
  final case class WorkingWithGcsData(requester: ActorRef, action: BuildWorkflowMetadataJsonAction, gcsResponse: Json) extends WorkingData

  implicit class EnhancedJson(val json: Json) extends AnyVal {
    def editFor(action: BuildWorkflowMetadataJsonAction): ErrorOr[Json] = action match {
      case _: GetLogs => JsonEditor.logs(json).validNel
      case _: WorkflowOutputs => JsonEditor.outputs(json).validNel
      case get: GetMetadataAction =>
        val intermediate = get.key match {
          case MetadataQuery(_, None, None, None, None, _) => json
          case MetadataQuery(_, None, Some(key), None, None, _) => JsonEditor.includeJson(json, NonEmptyList.of(key))
          case MetadataQuery(_, None, None, includeKeys, excludeKeys, _) => JsonEditor.includeExcludeJson(json, includeKeys, excludeKeys)
          // The following three don't bother with exclusions or inclusions since they are only used internally
          // and the existing client code is robust to superfluous data.
          case MetadataQuery(_, Some(_), None, None, None, _) => json
          case MetadataQuery(_, Some(_), Some(_), None, None, _) => json
          case MetadataQuery(_, Some(_), None, _, _, _) => json
          case wrong => throw new RuntimeException(s"Programmer Error: Invalid MetadataQuery: $wrong")
        }
        // For carbonited metadata, "expanded" subworkflows translates to not deleting subworkflows out of the root workflow that already
        // contains them. So `_.validNel` for expanded subworkflows and `JsonEditor.removeSubworkflowMetadata` for unexpanded subworkflows.
        val processSubworkflowMetadata: Json => ErrorOr[Json] = if (get.key.expandSubWorkflows) _.validNel else JsonEditor.replaceSubworkflowMetadataWithId
        processSubworkflowMetadata(intermediate)
      case other =>
        throw new RuntimeException(s"Programmer Error: Unexpected BuildWorkflowMetadataJsonAction message of type '${other.getClass.getSimpleName}': $other")
    }

    def updateLabels(labels: Map[WorkflowId, Map[String, String]]): ErrorOr[Json] = JsonEditor.updateLabels(json, labels)

    def toSpray: JsObject = json.printWith(Printer.spaces2).parseJson.asJsObject
  }

  implicit class EnhancedBuildWorkflowMetadataJsonAction(val action: BuildWorkflowMetadataJsonAction) extends AnyVal {
    def requiresLabels: Boolean = action match {
      case _: GetLogs | _: WorkflowOutputs => false
      case q: GetMetadataAction if q.key.key.isDefined && !q.key.key.contains("labels") => false
      case q: GetMetadataAction if q.key.excludeKeysOption.exists { _.toList.contains("labels") } => false
      case q: GetMetadataAction if q.key.includeKeysOption.exists { _.toList.contains("labels") } => true
      case q: GetMetadataAction if q.key.includeKeysOption.isDefined => false
      case _ => true
    }
  }

  implicit class EnhancedRequester(val requester: ActorRef) extends AnyVal {
    private def respondWith(action: BuildWorkflowMetadataJsonAction, sprayJson: ErrorOr[JsObject]): Unit = {
      val message = sprayJson match {
        case Valid(json) => SuccessfulMetadataJsonResponse(action, json)
        case Invalid(errors) => FailedMetadataJsonResponse(action, new Throwable(errors.toList.mkString("\n")))
      }
      requester ! message
    }

    def respondWith(action: BuildWorkflowMetadataJsonAction, gcsResponse: Json): Unit = {
      val responseJson = gcsResponse.editFor(action) map { _.toSpray }
      requester.respondWith(action, responseJson)
    }

    def respondWith(action: BuildWorkflowMetadataJsonAction, gcsResponse: Json, labels: Map[WorkflowId, Map[String, String]]): Unit = {
      val responseJson = for {
        edited <- gcsResponse.editFor(action)
        updated <- edited.updateLabels(labels)
      } yield updated.toSpray

      requester.respondWith(action, responseJson)
    }
  }
}
