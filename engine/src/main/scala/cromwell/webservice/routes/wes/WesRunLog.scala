package cromwell.webservice.routes.wes

import spray.json.JsObject

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
