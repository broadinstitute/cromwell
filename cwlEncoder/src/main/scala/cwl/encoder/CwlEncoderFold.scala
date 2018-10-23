package cwl.encoder

import cwl.CwlCodecs._
import cwl.{CommandLineTool, ExpressionTool, Workflow}
import io.circe.Json
import io.circe.syntax._
import shapeless.Poly1

object CwlEncoderFold extends Poly1 {
  implicit val workflow: Case.Aux[Workflow, Json] = at { _.asJson }
  implicit val commandLineTool: Case.Aux[CommandLineTool, Json] = at { _.asJson }
  implicit val expressionTool: Case.Aux[ExpressionTool, Json] = at { _.asJson }
}
