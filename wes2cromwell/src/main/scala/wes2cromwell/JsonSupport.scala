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
  implicit val workflowDescriptionFormat = jsonFormat2(WorkflowDescription)
  implicit val workflowLogEntryFormat = jsonFormat7(WorkflowLogEntry)
  implicit val workflowRequestFormat = jsonFormat7(WorkflowRequest)
  implicit val workflowLogFormat = jsonFormat6(WorkflowLog)
  implicit val workflowTypeVersionFormat = jsonFormat1(WorkflowTypeVersion)
  implicit val errorResponseFormat = jsonFormat2(ErrorResponse)
  implicit val wesResponseErrorFormat = jsonFormat2(WesResponseError)
  implicit val wesResponseCreateWorkflowIdFormat = jsonFormat1(WesResponseCreateWorkflowId)
  implicit val wesResponseDeleteWorkflowIdFormat = jsonFormat1(WesResponseDeleteWorkflowId)
  implicit val wesResponseStatusFormat = jsonFormat2(WesResponseStatus)
  implicit val wesResponseWorkflowListFormat = jsonFormat1(WesResponseWorkflowList)
}
