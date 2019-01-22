package cromwell.services.womtool.models

import io.circe.Decoder.Result
import io.circe.{Decoder, Encoder, HCursor, Json}
import wom.expression.WomExpression

object WomExpressionJsonSupport {
  implicit val womExpressionEncoder: Encoder[WomExpression] = new Encoder[WomExpression] {
    override def apply(a: WomExpression): Json = Json.fromString(a.sourceString)
  }

  implicit val womExpressionDecoder: Decoder[WomExpression] = new Decoder[WomExpression] {
    override def apply(c: HCursor): Result[WomExpression] = ???
  }
}
