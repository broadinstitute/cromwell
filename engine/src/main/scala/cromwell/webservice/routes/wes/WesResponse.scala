package cromwell.webservice.routes.wes

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import cromwell.webservice.routes.wes.WesState.WesState
import spray.json.{DefaultJsonProtocol, JsObject, JsonFormat, RootJsonFormat}

sealed trait WesResponse extends Product with Serializable
final case class WesErrorResponse(msg: String, status_code: Int) extends WesResponse
final case class WesRunId(run_id: String) extends WesResponse
final case class WesRunStatus(run_id: String, state: WesState) extends WesResponse

final case class WesResponseRunLog(run_id: String,
                                   request: WesRunRequest,
                                   state: WesState,
                                   run_log: Option[WesLog],
                                   task_logs: Option[List[WesLog]],
                                   outputs: Option[JsObject]) extends WesResponse

object WesResponseRunLog {
  def fromJson(json: JsObject): WesResponseRunLog = CromwellMetadata.fromJson(json).wesResponseRunLog
}

final case class WesLog(name: Option[String],
                        cmd: Option[Seq[String]],
                        start_time: Option[String],
                        end_time: Option[String],
                        stdout: Option[String],
                        stderr: Option[String],
                        exit_code: Option[Int]
                       )

final case class WesRunRequest(workflow_params: Option[JsObject],
                               workflow_type: Option[String],
                               workflow_type_version: Option[String],
                               tags: Option[JsObject],
                               workflow_engine_parameters: Option[JsObject],
                               workflow_url: Option[String]
                              )

final case class WesStatusInfoResponse(workflow_type_version: Map[String, Iterable[String]],
                                    supported_wes_versions: Iterable[String],
                                    supported_filesystem_protocols: Iterable[String],
                                    workflow_engine_versions: Map[String, String],
                                    default_workflow_engine_parameters: Iterable[DefaultWorkflowEngineParameter],
                                    system_state_counts: Map[WesState, Int],
                                    auth_instructions_url: String,
                                    contact_info_url: String,
                                    tags: Map[String, String]) extends WesResponse

object WesResponseJsonSupport extends SprayJsonSupport with DefaultJsonProtocol {
  import WesStateJsonSupport._
  import DefaultWorkflowEngineParameter.DefaultWorkflowEngineParameterFormat

  implicit val WesResponseErrorFormat = jsonFormat2(WesErrorResponse)
  implicit val WesResponseRunIdFormat = jsonFormat1(WesRunId)
  implicit val WesResponseStatusFormat = jsonFormat2(WesRunStatus)
  implicit val WesResponseStatusInfoFormat = jsonFormat9(WesStatusInfoResponse)
  implicit val WesLogFormat: JsonFormat[WesLog] = jsonFormat7(WesLog)
  implicit val WesRunRequestFormat: JsonFormat[WesRunRequest] = jsonFormat6(WesRunRequest)
  implicit val WesResponseRunLogFormat = jsonFormat6(WesResponseRunLog.apply)

  implicit object WesResponseFormat extends RootJsonFormat[WesResponse] {
    import spray.json._

    def write(r: WesResponse) = {
      r match {
        case r: WesRunId => r.toJson
        case s: WesRunStatus => s.toJson
        case e: WesErrorResponse => e.toJson
        case l: WesResponseRunLog => l.toJson
        case i: WesStatusInfoResponse => i.toJson
      }
    }

    def read(value: JsValue) = throw new UnsupportedOperationException("Reading WesResponse objects from JSON is not supported")
  }
}
