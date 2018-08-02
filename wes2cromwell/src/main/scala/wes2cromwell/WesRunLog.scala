package wes2cromwell

import spray.json.{JsObject, JsonFormat}

final case class WesLog(name: Option[String],
                        cmd: Option[Seq[String]],
                        start_time: Option[String],
                        end_time: Option[String],
                        stdout: Option[String],
                        stderr: Option[String],
                        exit_code: Option[Int]
                       )

final case class WesRunRequest(workflow_params: Option[JsObject],
                               workflow_type: String,
                               workflow_type_version: String,
                               tags: Option[JsObject],
                               workflow_engine_parameters: Option[JsObject],
                               workflow_url: Option[String]
                              )

final case class WesRunLog(run_id: String,
                           request: WesRunRequest,
                           state: WesState,
                           run_log: Option[WesLog],
                           task_logs: Option[List[WesLog]],
                           outputs: Option[JsObject]
                          )

object WesRunLog {
  def fromJson(json: String): WesRunLog = CromwellMetadata.fromJson(json).wesRunLog
}

object WorkflowLogJsonSupport {
  import WesStateJsonSupport._

  implicit val logFormat: JsonFormat[WesLog] = jsonFormat7(WesLog)
  implicit val runRequestFormat: JsonFormat[WesRunRequest] = jsonFormat6(WesRunRequest)
  implicit val runLogFormat: JsonFormat[WesRunLog] = jsonFormat6(WesRunLog.apply)
}
