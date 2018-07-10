package wes2cromwell

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import spray.json.{ DefaultJsonProtocol, JsString, JsValue, RootJsonFormat }

trait JsonSupport extends SprayJsonSupport {
  // import the default encoders for primitive types (Int, String, Lists etc)
  import DefaultJsonProtocol._

  // Json marshall/unmarshall for WorkflowState
  // The order of implicits seems to be important; they have to come in reference order.
  // TODO: tried and failed to get this code to live with WorkflowState.scala. Never figured out why.
  implicit object WorkflowStateFormat extends RootJsonFormat[WorkflowState] {
    def write(obj: WorkflowState): JsValue = JsString(obj.toString)

    // TODO: not sure this works and nothing in the code tests it right now
    def read(json: JsValue): WorkflowState = {
      json match {
        case x: JsString => WorkflowState.withName(x.value)
        case _ => throw new IllegalArgumentException("Expected a string value for state")
      }
    }
  }
  // WES structures
  implicit val workflowDescriptionFormat: RootJsonFormat[WorkflowDescription] = jsonFormat2(WorkflowDescription)
  implicit val workflowLogEntryFormat: RootJsonFormat[WorkflowLogEntry] = jsonFormat7(WorkflowLogEntry)
  implicit val workflowRequestFormat: RootJsonFormat[WesSubmission] = jsonFormat7(WesSubmission)
  implicit val workflowLogFormat: RootJsonFormat[WorkflowLog] = jsonFormat6(WorkflowLog)
  implicit val workflowTypeVersionFormat: RootJsonFormat[WorkflowTypeVersion] = jsonFormat1(WorkflowTypeVersion)
  implicit val errorResponseFormat: RootJsonFormat[ErrorResponse] = jsonFormat2(ErrorResponse)
  implicit val wesResponseErrorFormat: RootJsonFormat[WesResponseError] = jsonFormat2(WesResponseError)
  implicit val wesResponseCreateWorkflowIdFormat: RootJsonFormat[WesResponseCreateWorkflowId] = jsonFormat1(WesResponseCreateWorkflowId)
  implicit val wesResponseDeleteWorkflowIdFormat: RootJsonFormat[WesResponseDeleteWorkflowId] = jsonFormat1(WesResponseDeleteWorkflowId)
  implicit val wesResponseStatusFormat: RootJsonFormat[WesResponseStatus] = jsonFormat2(WesResponseStatus)
  implicit val wesResponseWorkflowListFormat: RootJsonFormat[WesResponseWorkflowList] = jsonFormat1(WesResponseWorkflowList)
  implicit val WesResponseWorkflowMetadataFormat: RootJsonFormat[WesResponseWorkflowMetadata] = jsonFormat1(WesResponseWorkflowMetadata)
}
