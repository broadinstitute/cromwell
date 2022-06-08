package cromwell.webservice.routes.wes

import cromwell.webservice.routes.wes.WesState.WesState
import spray.json.{JsObject, JsonParser}


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

final case class CromwellCallsMetadata(shardIndex: Option[Int],
                                       commandLine: Option[String],
                                       returnCode: Option[Int],
                                       start: Option[String],
                                       end: Option[String],
                                       stdout: Option[String],
                                       stderr: Option[String]
                                      )

object WesRunLog {
  def fromJson(json: String): WesRunLog = CromwellMetadata.fromJson(json).wesRunLog
}
