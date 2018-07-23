package wes2cromwell

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}


sealed trait NewWesRunState extends Product with Serializable
case object UNKNOWN extends NewWesRunState
case object QUEUED extends NewWesRunState
case object INITIALIZING extends NewWesRunState
case object RUNNING extends NewWesRunState
case object PAUSED extends NewWesRunState
case object COMPLETE extends NewWesRunState
case object EXECUTOR_ERROR extends NewWesRunState
case object SYSTEM_ERROR extends NewWesRunState
case object CANCELED extends NewWesRunState

object NewWesRunState {
  def fromCromwellState(cromwellState: String): NewWesRunState = {
    cromwellState match {
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

object NewWesRunStateJsonSupport extends SprayJsonSupport with DefaultJsonProtocol {
  implicit object RunStateFormat extends RootJsonFormat[NewWesRunState] {
    def write(obj: NewWesRunState): JsValue = JsString(obj.toString)
    def read(json: JsValue): NewWesRunState = throw new UnsupportedOperationException("Reading WES RunState unsupported")
  }
}

sealed trait WesRunState

// FIXME
object WesRunState {
  lazy val RunStateValues = Seq(UNKNOWN, QUEUED, INITIALIZING, RUNNING, PAUSED, COMPLETE, EXECUTOR_ERROR, SYSTEM_ERROR, CANCELED)

  implicit object RunStateFormat extends RootJsonFormat[WesRunState] {
    def write(obj: WesRunState): JsValue = JsString(obj.toString)
    def read(json: JsValue): WesRunState = throw new UnsupportedOperationException("Reading WES RunState unsupported")
  }

  case object UNKNOWN extends WesRunState {
    override val toString: String = "UNKNOWN"
  }
  case object QUEUED extends WesRunState {
    override val toString: String = "QUEUED"
  }
  case object INITIALIZING extends WesRunState {
    override val toString: String = "INITIALIZING"
  }
  case object RUNNING extends WesRunState {
    override val toString = "RUNNING"
  }
  case object PAUSED extends WesRunState {
    override val toString = "PAUSED"
  }
  case object COMPLETE extends WesRunState {
    override val toString = "COMPLETE"
  }
  case object EXECUTOR_ERROR extends WesRunState {
    override val toString = "EXECUTOR_ERROR"
  }
  case object SYSTEM_ERROR extends WesRunState {
    override val toString = "SYSTEM_ERROR"
  }
  case object CANCELED extends WesRunState {
    override val toString = "CANCELED"
  }
}

