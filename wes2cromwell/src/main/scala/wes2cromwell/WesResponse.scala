package wes2cromwell

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import spray.json.{DefaultJsonProtocol, RootJsonFormat}

// The WesResponse provides a trait for all possible responses to requests to the WES REST API
sealed trait WesResponse extends Product with Serializable

final case class WesErrorResponse(msg: String, status_code: Int) extends WesResponse
final case class WesRunId(run_id: String) extends WesResponse
final case class WesRunStatus(run_id: String, state: NewWesRunState) extends WesResponse
final case class WesResponseRunList(runs: List[WesRunStatus]) extends WesResponse
final case class WesResponseWorkflowMetadata(workflowLog: WesRunLog) extends WesResponse

object WesResponseJsonSupport extends SprayJsonSupport with DefaultJsonProtocol {
  import WorkflowLogJsonSupport._
  import NewWesRunStateJsonSupport._

  implicit val WesResponseErrorFormat = jsonFormat2(WesErrorResponse)
  implicit val WesResponseRunIdFormat = jsonFormat1(WesRunId)
  implicit val WesResponseStatusFormat = jsonFormat2(WesRunStatus)
  implicit val WesResponseRunListFormat = jsonFormat1(WesResponseRunList)
  implicit val WesResponseRunMetadataFormat = jsonFormat1(WesResponseWorkflowMetadata)

  implicit object WesResponseFormat extends RootJsonFormat[WesResponse] {
    import spray.json._

    def write(r: WesResponse) = {
      r match {
        case r: WesRunId => r.toJson
        case s: WesRunStatus => s.toJson
        case l: WesResponseRunList => l.toJson
        case e: WesErrorResponse => e.toJson
        case m: WesResponseWorkflowMetadata => m.toJson
      }
    }

    def read(value: JsValue) = throw new UnsupportedOperationException("Reading WesResponse objects from JSON is not supported")
  }
}

