package cromwell.webservice.routes.wes

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import cromwell.core._
import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}

object WesState {
  sealed trait WesState extends Product with Serializable
  case object UNKNOWN extends WesState
  case object QUEUED extends WesState
  case object INITIALIZING extends WesState
  case object RUNNING extends WesState
  case object PAUSED extends WesState
  case object COMPLETE extends WesState
  case object EXECUTOR_ERROR extends WesState
  case object SYSTEM_ERROR extends WesState
  case object CANCELED extends WesState
  case object CANCELING extends WesState // Before any shorebirds worry about the spelling, it really is CANCELING

  def fromCromwellStatus(cromwellStatus: WorkflowState): WesState = {
    cromwellStatus match {
      case WorkflowOnHold => PAUSED
      case WorkflowSubmitted => QUEUED
      case WorkflowRunning => RUNNING
      case WorkflowAborting => CANCELING
      case WorkflowAborted => CANCELED
      case WorkflowSucceeded => COMPLETE
      case WorkflowFailed => EXECUTOR_ERROR
      case _ => UNKNOWN
    }
  }
}

object WesStateJsonSupport extends SprayJsonSupport with DefaultJsonProtocol {
  import WesState._

  implicit object WesStateFormat extends RootJsonFormat[WesState] {
    def write(obj: WesState): JsValue = JsString(obj.toString)

    def read(json: JsValue): WesState = json match {
      case JsString(string) if string == "UNKNOWN" => UNKNOWN
      case JsString(string) if string == "QUEUED" => QUEUED
      case JsString(string) if string == "INITIALIZING" => INITIALIZING
      case JsString(string) if string == "RUNNING" => RUNNING
      case JsString(string) if string == "PAUSED" => PAUSED
      case JsString(string) if string == "COMPLETE" => COMPLETE
      case JsString(string) if string == "EXECUTOR_ERROR" => EXECUTOR_ERROR
      case JsString(string) if string == "SYSTEM_ERROR" => SYSTEM_ERROR
      case JsString(string) if string == "CANCELED" => CANCELED
      case JsString(string) if string == "CANCELING" => CANCELING
      case other => throw new UnsupportedOperationException(s"Cannot deserialize $other into a WesState")
    }
  }
}
