package wes2cromwell

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import spray.json.{DefaultJsonProtocol, JsString, JsValue, JsonFormat, RootJsonFormat}

// FIXME: REMOVE
trait JsonSupport extends SprayJsonSupport {
  // import the default encoders for primitive types (Int, String, Lists etc)
  import DefaultJsonProtocol._

  // Json marshall/unmarshall for WorkflowState
  // The order of implicits seems to be important; they have to come in reference order.
  // TODO: tried and failed to get this code to live with WorkflowState.scala. Never figured out why.
  implicit object WorkflowStateFormat extends RootJsonFormat[WesState] {
    def write(obj: WesState): JsValue = JsString(obj.toString)
    def read(json: JsValue): WesState = throw new UnsupportedOperationException("Reading WES RunState unsupported")
  }
}
