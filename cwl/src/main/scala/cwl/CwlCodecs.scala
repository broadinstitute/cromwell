package cwl

import io.circe._
import io.circe.parser._
import io.circe.generic.auto._
import eu.timepit.refined.string._
import io.circe.refined._
import io.circe.literal._
import common.Checked

object CwlCodecs {

  import cwl.decoder._
  implicit val cwlTypeDecoder = Decoder.enumDecoder(CwlType)
  implicit val cwlVersionDecoder = Decoder.enumDecoder(CwlVersion)
  implicit val scatterMethodDecoder = Decoder.enumDecoder(ScatterMethod)
  implicit val linkMergeMethodDecoder = Decoder.enumDecoder(LinkMergeMethod)

  //According to automatic derivation, these instances should not be required.  But
  //removing these breaks decodeCwl, so...
  implicit val wfD = implicitly[Decoder[Workflow]]
  implicit val cltD = implicitly[Decoder[CommandLineTool]]
  implicit val etD = implicitly[Decoder[ExpressionTool]]

  def decodeCwl(in: String): Checked[CwlFile] = decodeAccumulating[CwlFile](in).leftMap(_.map(_.getMessage)).toEither
}
