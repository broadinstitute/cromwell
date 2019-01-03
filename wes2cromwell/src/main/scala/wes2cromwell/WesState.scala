package wes2cromwell

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}

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

object WesState {
  def fromCromwellStatus(cromwellStatus: String): WesState = {
    cromwellStatus match {
      case "On Hold" => PAUSED
      case "Submitted" => QUEUED
      case "Running" => RUNNING
      case "Aborting" => CANCELED
      case "Aborted" => CANCELED
      case "Succeeded" => COMPLETE
      case "Failed" => EXECUTOR_ERROR
      case _ => UNKNOWN
    }
  }
}

object WesStateJsonSupport extends SprayJsonSupport with DefaultJsonProtocol {
  implicit object WesStateFormat extends RootJsonFormat[WesState] {
    def write(obj: WesState): JsValue = JsString(obj.toString)

    def read(json: JsValue): WesState = throw new UnsupportedOperationException("Reading WesState unsupported")
  }
}
