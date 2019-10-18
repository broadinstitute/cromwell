package cromwell.services.metadata.hybridcarbonite

import akka.actor.{ActorRef, LoggingFSM, Props, Status}
import akka.pattern.pipe
import cromwell.core.io.{AsyncIo, DefaultIoCommandBuilder}
import cromwell.core.{RootWorkflowId, WorkflowId}
import cromwell.services.metadata.MetadataService.{GetRootAndSubworkflowLabels, RootAndSubworkflowLabelsLookupFailed, RootAndSubworkflowLabelsLookupResponse}
import cromwell.services.metadata.hybridcarbonite.CarbonitedMetadataThawingActor._
import cromwell.util.JsonEditor
import io.circe.parser._
import io.circe.{Json, Printer}

import scala.concurrent.{ExecutionContext, Future}

final class CarbonitedMetadataThawingActor(carboniterConfig: HybridCarboniteConfig, serviceRegistry: ActorRef, ioActor: ActorRef) extends LoggingFSM[CarbonitedMetadataThawingState, CarbonitedMetadataThawingData] {

  implicit val ec: ExecutionContext = context.dispatcher

  val asyncIo = new AsyncIo(ioActor, DefaultIoCommandBuilder)

  startWith(CarbonitedMetadataThawingActorPendingState, CarbonitedMetadataThawingActorPendingData)

  when(CarbonitedMetadataThawingActorPendingState) {
    case Event(ThawCarbonitedMetadata(rootWorkflowId), CarbonitedMetadataThawingActorPendingData) =>

      val response = for {
        rawResponseFromIoActor <- asyncIo.contentAsStringAsync(carboniterConfig.makePath(rootWorkflowId), maxBytes = Option(30 * 1024 * 1024), failOnOverflow = true)
        parsedResponse <- Future.fromTry(parse(rawResponseFromIoActor).toTry.map(GcsMetadataResponse.apply))
      } yield parsedResponse

      response pipeTo self
      serviceRegistry ! GetRootAndSubworkflowLabels(rootWorkflowId)

      goto(CarbonitedMetadataThawingActorRunningState) using CarbonitedMetadataThawingActorDataNothingYet(sender())
  }

  when(CarbonitedMetadataThawingActorRunningState) {

    // Successful Responses:
    case Event(GcsMetadataResponse(response), CarbonitedMetadataThawingActorDataNothingYet(requester)) =>
      stay() using CarbonitedMetadataThawingActorDataWithGcsResponse(requester, response)

    case Event(GcsMetadataResponse(response), CarbonitedMetadataThawingActorDataWithLabels(requester, labels)) =>
      requester ! ThawCarboniteSucceeded(mergeResponses(response, labels))
      stop()

    case Event(RootAndSubworkflowLabelsLookupResponse(_, labels), CarbonitedMetadataThawingActorDataNothingYet(requester)) =>
      stay() using CarbonitedMetadataThawingActorDataWithLabels(requester, labels)

    case Event(RootAndSubworkflowLabelsLookupResponse(_, labels), CarbonitedMetadataThawingActorDataWithGcsResponse(requester, gcsResponse)) =>
      requester ! ThawCarboniteSucceeded(mergeResponses(gcsResponse, labels))
      stop()

    // Failures:
    case Event(Status.Failure(failure), data: CarbonitedMetadataThawingActorWorkingData) =>
      data.requester ! ThawCarboniteFailed(failure)
      stop()
    case Event(RootAndSubworkflowLabelsLookupFailed(_, failure), data: CarbonitedMetadataThawingActorWorkingData) =>
      data.requester ! ThawCarboniteFailed(failure)
      stop()
  }

  whenUnhandled {
    case other =>
      log.error(s"Programmer Error: Unexpected message to ${self.path.name}: $other")
      stay()
  }

  def mergeResponses(metadata: Json, labels: Map[WorkflowId, Map[String, String]]): String = JsonEditor.updateLabels(metadata, labels).printWith(Printer.spaces2)

}

object CarbonitedMetadataThawingActor {

  def props(carboniterConfig: HybridCarboniteConfig, serviceRegistry: ActorRef, ioActor: ActorRef): Props =
    Props(new CarbonitedMetadataThawingActor(carboniterConfig, serviceRegistry, ioActor))

  sealed trait ThawCarboniteMessage

  final case class ThawCarbonitedMetadata(workflowId: RootWorkflowId) extends ThawCarboniteMessage
  final case class ThawCarboniteSucceeded(value: String) extends ThawCarboniteMessage
  final case class ThawCarboniteFailed(reason: Throwable) extends ThawCarboniteMessage

  final case class GcsMetadataResponse(response: Json)

  sealed trait CarbonitedMetadataThawingState
  case object CarbonitedMetadataThawingActorPendingState extends CarbonitedMetadataThawingState
  case object CarbonitedMetadataThawingActorRunningState extends CarbonitedMetadataThawingState

  sealed trait CarbonitedMetadataThawingData
  case object CarbonitedMetadataThawingActorPendingData extends CarbonitedMetadataThawingData
  sealed trait CarbonitedMetadataThawingActorWorkingData extends CarbonitedMetadataThawingData { def requester: ActorRef }
  final case class CarbonitedMetadataThawingActorDataNothingYet(requester: ActorRef) extends CarbonitedMetadataThawingActorWorkingData
  final case class CarbonitedMetadataThawingActorDataWithGcsResponse(requester: ActorRef, gcsResponse: Json) extends CarbonitedMetadataThawingActorWorkingData
  final case class CarbonitedMetadataThawingActorDataWithLabels(requester: ActorRef, labels: Map[WorkflowId, Map[String, String]]) extends CarbonitedMetadataThawingActorWorkingData

}
