package cwl

import cats.data.NonEmptyList
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

  def decodeCwl(in: String): Checked[CwlFile] = {
    decodeAccumulating[CwlFile](in).leftMap((original: NonEmptyList[Error]) => {

      def errorToStrings : Error => NonEmptyList[String] = e => NonEmptyList(e.getMessage, e.getStackTrace.map(_.toString).toList)

      val originalErrorMessageAsStrings: NonEmptyList[String] = original.flatMap(errorToStrings)
      //we know something is wrong, but we'd like to understand it better
      def betterError: String => NonEmptyList[String] = s => (s match {
        case "Workflow" => decode[Workflow](in)
        case "CommandLineTool" => decode[CommandLineTool](in)
        case "ExpressionTool" => decode[ExpressionTool](in)
        case _ => throw new Exception(s"saw an unknown CWL type: $s")
      }).fold(errorToStrings,
        _ => originalErrorMessageAsStrings
      )

      (for {
        raw <- parse(in).toOption
        clazz <- (raw \\ "class").head.asString
        errors = betterError(clazz)
      } yield errors).getOrElse(originalErrorMessageAsStrings)
    }
    ).toEither
  }
}
