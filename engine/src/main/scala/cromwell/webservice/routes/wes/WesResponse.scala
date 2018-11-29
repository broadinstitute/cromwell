package cromwell.webservice.routes.wes

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import cromwell.webservice.routes.wes.WesState.WesState
import spray.json.{DefaultJsonProtocol, RootJsonFormat}

sealed trait WesResponse extends Product with Serializable
final case class WesErrorResponse(msg: String, status_code: Int) extends WesResponse
final case class WesRunId(run_id: String) extends WesResponse
final case class WesRunStatus(run_id: String, state: WesState) extends WesResponse

object WesResponseJsonSupport extends SprayJsonSupport with DefaultJsonProtocol {
  import WesStateJsonSupport._

  implicit val WesResponseErrorFormat = jsonFormat2(WesErrorResponse)
  implicit val WesResponseRunIdFormat = jsonFormat1(WesRunId)
  implicit val WesResponseStatusFormat = jsonFormat2(WesRunStatus)

  implicit object WesResponseFormat extends RootJsonFormat[WesResponse] {
    import spray.json._

    def write(r: WesResponse) = {
      r match {
        case r: WesRunId => r.toJson
        case s: WesRunStatus => s.toJson
        case e: WesErrorResponse => e.toJson
      }
    }

    def read(value: JsValue) = throw new UnsupportedOperationException("Reading WesResponse objects from JSON is not supported")
  }
}
