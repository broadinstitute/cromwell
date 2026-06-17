package cromwell.engine.workflow.lifecycle.finalization

import cromwell.util.JsonFormatting.WomValueJsonFormatter.WomValueJsonFormat
import spray.json.DefaultJsonProtocol
import spray.json.RootJsonFormat
import wom.values.WomValue

final case class CallbackMessage(workflowId: String,
                                 state: String,
                                 outputs: Map[String, WomValue],
                                 failures: List[String]
)

object WorkflowCallbackJsonSupport extends DefaultJsonProtocol {
  implicit val callbackMessageFormat: RootJsonFormat[CallbackMessage] = jsonFormat4(CallbackMessage)
}
