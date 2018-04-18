package cwl

import cats.data.NonEmptyList
import io.circe._
import io.circe.parser._
import io.circe.generic.auto._
import eu.timepit.refined.string._
import io.circe.refined._
import io.circe.literal._
import common.Checked
import common.validation.Checked._
import cats.syntax.either._
import cats.syntax.show._
import cwl.CwlType.CwlType
import cwl.CwlVersion.CwlVersion
import cwl.LinkMergeMethod.LinkMergeMethod
import cwl.ScatterMethod.ScatterMethod
import io.circe.Json._
import io.circe.DecodingFailure._
import shapeless.Coproduct

object CwlCodecs {

  import cwl.decoder._
  implicit val cwlTypeDecoder        : Decoder[CwlType]         = Decoder.enumDecoder(CwlType)
  implicit val cwlVersionDecoder     : Decoder[CwlVersion]      = Decoder.enumDecoder(CwlVersion)
  implicit val scatterMethodDecoder  : Decoder[ScatterMethod]   = Decoder.enumDecoder(ScatterMethod)
  implicit val linkMergeMethodDecoder: Decoder[LinkMergeMethod] = Decoder.enumDecoder(LinkMergeMethod)

  //According to automatic derivation, these instances should not be required.  But
  //removing these breaks decodeCwl, so...
  implicit lazy val wfD : Decoder[Workflow] = implicitly[Decoder[Workflow]]
  implicit lazy val cltD: Decoder[CommandLineTool] = implicitly[Decoder[CommandLineTool]]
  implicit lazy val etD : Decoder[ExpressionTool] = implicitly[Decoder[ExpressionTool]]

  def decodeCwl(json: Json): Checked[Cwl] = {
    findClass(json) match {
      case Some("Workflow") => decodeWithErrorStringsJson[Workflow](json).map(Coproduct[Cwl].apply(_))
      case Some("CommandLineTool") => decodeWithErrorStringsJson[CommandLineTool](json).map(Coproduct[Cwl].apply(_))
      case Some("ExpressionTool") => decodeWithErrorStringsJson[ExpressionTool](json).map(Coproduct[Cwl].apply(_))
      case Some(other) => s"Class field was declared incorrectly: $other is not one of Workflow, CommandLineTool, or ExpressionTool! as seen in ${json.show}".invalidNelCheck
      case None => s"Class field was omitted in ${json.show}".invalidNelCheck
    }
  }

  private def decodeWithErrorStringsJson[A](in: Json)(implicit d: Decoder[A]): Checked[A] =
    in.as[A].leftMap(_.show).leftMap(NonEmptyList.one)

  private def findClass(json: Json): Option[String] =
    for {
      obj <- json.asObject
      map = obj.toMap
      classObj <- map.get("class")
      classString <- classObj.asString
    } yield classString
}
