package cromwell.services.womtool.models

import spray.json.{JsonWriter, JsString, JsValue}
import wom.expression.WomExpression

object WomExpressionJsonSupport {
  implicit val womExpressionEncoder: JsonWriter[WomExpression] = new JsonWriter[WomExpression] {
    override def write(a: WomExpression): JsValue = JsString(a.sourceString)
  }
}
