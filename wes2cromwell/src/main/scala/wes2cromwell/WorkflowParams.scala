package wes2cromwell

import java.util.Base64

import spray.json.{ DefaultJsonProtocol, JsonFormat, JsonParser }

/*
 * The WorkflowParams class is in the structure we expect to see in the JSON coming in the
 * workflow_params field in the POST /workflows REST API. The object contains a factory method
 * to generate the class from a JSON string
 */
final case class WorkflowParams(
    workflowOnHold: Option[Boolean],
    workflowInputs: List[String],
    workflowOptions: Option[String],
    workflowDependencies: Option[String]
) {

  def dependenciesZip(): Option[String] = {
    workflowDependencies match {
      case Some(x) => Some(Base64.getDecoder.decode(x).toString)
      case None => None
    }
  }
}

object WorkflowParams {
  // WES-Opaque structure
  import DefaultJsonProtocol._
  implicit val workflowParamsFormat: JsonFormat[WorkflowParams] = jsonFormat4(WorkflowParams.apply)

  def toWorkflowParams(json: String): WorkflowParams = {
    val jsonAst = JsonParser(json)
    jsonAst.convertTo[WorkflowParams]
  }
}
