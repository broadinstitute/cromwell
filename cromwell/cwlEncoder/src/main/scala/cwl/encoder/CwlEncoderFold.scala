package cwl.encoder

import cwl.{CommandLineTool, CwlType, CwlVersion, ExpressionTool, LinkMergeMethod, ScatterMethod, Workflow}
import io.circe.{Encoder, Json}
import shapeless.Poly1
import io.circe._
import io.circe.shapes._
import io.circe.generic.auto._
import io.circe.refined._
import io.circe.literal._
import shapeless.Poly1
import cwl._
import io.circe.syntax._

object CwlEncoderFold extends Poly1 {

  implicit val cwlTypeEncoder = Encoder.enumEncoder(CwlType)
  implicit val cwlVersionEncoder = Encoder.enumEncoder(CwlVersion)
  implicit val scatterMethodEncoder = Encoder.enumEncoder(ScatterMethod)
  implicit val linkMergeMethodEncoder = Encoder.enumEncoder(LinkMergeMethod)

  implicit def enumerationEncoder[V <: Enumeration#Value]: Encoder[V] = Encoder.instance((value: V) => Json.fromString(value.toString))

  implicit def workflow = at[Workflow]{ _.asJson }

  implicit def commandLineTool = at[CommandLineTool]{ _.asJson }

  implicit def expressionTool = at[ExpressionTool]{ _.asJson }

}
