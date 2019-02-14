package cromwell.services.womtool.models

import io.circe.{Encoder, Json}
import wom.expression.WomExpression

object WomExpressionJsonSupport {
  implicit val womExpressionEncoder: Encoder[WomExpression] = new Encoder[WomExpression] {
    override def apply(a: WomExpression): Json = Json.fromString(a.sourceString)
  }
}
