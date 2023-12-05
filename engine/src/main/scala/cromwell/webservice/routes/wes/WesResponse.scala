package cromwell.webservice.routes.wes

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import cromwell.webservice.routes.wes.WesState.WesState
import spray.json.{DefaultJsonProtocol, JsObject, RootJsonFormat}

sealed trait WesResponse extends Product with Serializable
final case class WesErrorResponse(msg: String, status_code: Int) extends WesResponse
final case class WesRunId(run_id: String) extends WesResponse
final case class WesRunStatus(run_id: String, state: WesState) extends WesResponse
final case class WesResponseRunList(runs: List[WesRunStatus]) extends WesResponse
final case class WesResponseWorkflowMetadata(workflowLog: WesRunLog) extends WesResponse

final case class WesRunLog(run_id: String,
                           request: WesRunRequest,
                           state: WesState,
                           run_log: Option[WesLog],
                           task_logs: Option[List[WesLog]],
                           outputs: Option[JsObject]
) extends WesResponse

object WesRunLog {
  def fromJson(json: String): WesRunLog = CromwellMetadata.fromJson(json).wesRunLog
}

final case class WesStatusInfoResponse(workflow_type_version: Map[String, Iterable[String]],
                                       supported_wes_versions: Iterable[String],
                                       supported_filesystem_protocols: Iterable[String],
                                       workflow_engine_versions: Map[String, String],
                                       default_workflow_engine_parameters: Iterable[DefaultWorkflowEngineParameter],
                                       system_state_counts: Map[WesState, Int],
                                       auth_instructions_url: String,
                                       contact_info_url: String,
                                       tags: Map[String, String]
) extends WesResponse

object WesResponseJsonSupport extends SprayJsonSupport with DefaultJsonProtocol {
  import WesStateJsonSupport._
  import DefaultWorkflowEngineParameter.DefaultWorkflowEngineParameterFormat

  implicit val WesResponseErrorFormat = jsonFormat2(WesErrorResponse)
  implicit val WesResponseRunIdFormat = jsonFormat1(WesRunId)
  implicit val WesResponseStatusFormat = jsonFormat2(WesRunStatus)
  implicit val WesResponseStatusInfoFormat = jsonFormat9(WesStatusInfoResponse)
  implicit val WesResponseRunListFormat = jsonFormat1(WesResponseRunList)
  implicit val WesRunRequestFormat = jsonFormat6(WesRunRequest)
  implicit val WesLogFormat = jsonFormat7(WesLog)
  implicit val WesRunLogFormat = jsonFormat6(WesRunLog.apply)
  implicit val WesResponseWorkflowMetadataFormat = jsonFormat1(WesResponseWorkflowMetadata)

  implicit object WesResponseFormat extends RootJsonFormat[WesResponse] {
    import spray.json._

    def write(r: WesResponse) =
      r match {
        case r: WesRunId => r.toJson
        case s: WesRunStatus => s.toJson
        case e: WesErrorResponse => e.toJson
        case i: WesStatusInfoResponse => i.toJson
        case l: WesResponseRunList => l.toJson
        case m: WesResponseWorkflowMetadata => m.toJson
        case w: WesRunLog => w.toJson
      }

    def read(value: JsValue) = throw new UnsupportedOperationException(
      "Reading WesResponse objects from JSON is not supported"
    )
  }
}
