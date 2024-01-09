package cromwell.engine.workflow.lifecycle.finalization

import cromwell.util.JsonFormatting.WomValueJsonFormatter.WomValueJsonFormat
import spray.json.DefaultJsonProtocol
import wom.values.WomValue

final case class CallbackMessage(workflowId: String,
                                 state: String,
                                 outputs: Map[String, WomValue],
                                 failures: List[String]
)

object WorkflowCallbackJsonSupport extends DefaultJsonProtocol {
  implicit val callbackMessageFormat = jsonFormat4(CallbackMessage)
}
