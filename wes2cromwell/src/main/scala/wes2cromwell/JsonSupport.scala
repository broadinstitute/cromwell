package wes2cromwell

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import spray.json.{DefaultJsonProtocol, JsString, JsValue, JsonFormat, RootJsonFormat}

trait JsonSupport extends SprayJsonSupport {
  // import the default encoders for primitive types (Int, String, Lists etc)
  import DefaultJsonProtocol._

  // Json marshall/unmarshall for WorkflowState
  // The order of implicits seems to be important; they have to come in reference order.
  // TODO: tried and failed to get this code to live with WorkflowState.scala. Never figured out why.
  implicit object WorkflowStateFormat extends RootJsonFormat[WesRunState] {
    def write(obj: WesRunState): JsValue = JsString(obj.toString)
    def read(json: JsValue): WesRunState = throw new UnsupportedOperationException("Reading WES RunState unsupported")
  }


  // WES structures
  implicit val workflowDescriptionFormat: JsonFormat[WorkflowDescription] = jsonFormat2(WorkflowDescription)
  implicit val workflowTypeVersionFormat: JsonFormat[WorkflowTypeVersion] = jsonFormat1(WorkflowTypeVersion)
  implicit val errorResponseFormat: JsonFormat[ErrorResponse] = jsonFormat2(ErrorResponse)
}
