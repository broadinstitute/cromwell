package cromwell.webservice.routes.wes

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import cromwell.core._
import spray.json.{DefaultJsonProtocol, JsObject, JsString, JsValue, RootJsonFormat}

object WesState {
  sealed trait WesState extends Product with Serializable { val name: String }
  case object Unknown extends WesState { override val name = "UNKNOWN"}
  case object Queued extends WesState { override val name = "QUEUED"}
  case object Initializing extends WesState { override val name = "INITIALIZING"}
  case object Running extends WesState { override val name = "RUNNING"}
  case object Paused extends WesState { override val name = "PAUSED"}
  case object Complete extends WesState { override val name = "COMPLETE"}
  case object ExecutorError extends WesState { override val name = "EXECUTOR_ERROR"}
  case object SystemError extends WesState { override val name = "SYSTEM_ERROR"}
  case object Canceled extends WesState { override val name = "CANCELED"}
  case object Canceling extends WesState { override val name = "CANCELING"}

  def fromCromwellStatus(cromwellStatus: WorkflowState): WesState = {
    cromwellStatus match {
      case WorkflowOnHold => Paused
      case WorkflowSubmitted => Queued
      case WorkflowRunning => Running
      case WorkflowAborting => Canceling
      case WorkflowAborted => Canceled
      case WorkflowSucceeded => Complete
      case WorkflowFailed => ExecutorError
      case _ => Unknown
    }
  }

  def fromCromwellStatusJson(jsonResponse: JsObject): WesState = {

    val statusString = jsonResponse.fields.get("status").collect {
      case str: JsString => str.value
    }.getOrElse(throw new IllegalArgumentException(s"Could not coerce Cromwell status response ${jsonResponse.compactPrint} into a valid WES status"))

    fromCromwellStatus(WorkflowState.withName(statusString))
  }

  def fromString(status: String): WesState = {
    status match {
      case Unknown.name => Unknown
      case Queued.name => Queued
      case Initializing.name => Initializing
      case Running.name => Running
      case Paused.name => Paused
      case Complete.name => Complete
      case ExecutorError.name => ExecutorError
      case SystemError.name => SystemError
      case Canceled.name => Canceled
      case Canceling.name => Canceling
      case doh => throw new IllegalArgumentException(s"Invalid status attempting to be coerced to WesState: $doh")
    }
  }
}

object WesStateJsonSupport extends SprayJsonSupport with DefaultJsonProtocol {
  import WesState._

  implicit object WesStateFormat extends RootJsonFormat[WesState] {
    def write(obj: WesState): JsValue = JsString(obj.name)

    def read(json: JsValue): WesState =
      json match {
      case JsString(string) => WesState.fromString(string)
      case other => throw new UnsupportedOperationException(s"Cannot deserialize $other into a WesState")
    }
  }
}
