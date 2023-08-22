package cromwell.engine.workflow.lifecycle.finalization

import cromwell.util.JsonFormatting.WomValueJsonFormatter.WomValueJsonFormat
import spray.json.DefaultJsonProtocol
import wom.values.WomValue

final case class CallbackMessage(workflowId: String, state: String, outputs: Map[String, WomValue])

object WorkflowCallbackJsonSupport extends DefaultJsonProtocol {
  implicit val callbackMessageFormat = jsonFormat3(CallbackMessage)
}
