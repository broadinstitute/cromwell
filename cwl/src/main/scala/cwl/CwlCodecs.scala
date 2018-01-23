package cwl

import cats.data.NonEmptyList
import io.circe._
import io.circe.parser._
import io.circe.generic.auto._
import eu.timepit.refined.string._
import io.circe.refined._
import io.circe.literal._
import common.Checked
import common.validation.ErrorOr.ErrorOr
import cats.syntax.either._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.syntax.show._
import cats.instances.vector._
import io.circe.Error._
import io.circe.ParsingFailure._

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

  def errorToStrings : NonEmptyList[Error] => NonEmptyList[String] = _.map(_.show)

  def decodeWithErrorStrings[A](in: String)(implicit d: Decoder[A]): ErrorOr[A] = decodeAccumulating[A](in).leftMap(errorToStrings)

  /**
    * Attempt to parse as an array of CWL or a single one.
    *
    * If it fails, attempt to parse each CWL individually and return a superset of those failures.
    */
  def decodeCwl(in: String): Checked[CwlFile] =
    decodeWithErrorStrings[CwlFile](in).leftMap(decodePieces(in)).toEither

  /**
    * Drop down to a lower level and try to parse each CWL individually, returning a superset of those failures.
    */
  def decodePieces(in: String)(original: NonEmptyList[String]):NonEmptyList[String]  =
    (for {
      raw <- parse(in).leftMap(_.show).leftMap(NonEmptyList.one)
      _ <-
        if (raw.isArray)
          raw.asArray.toVector.flatten.traverse(decodeCwl).toEither
        else
          decodeCwl(raw).toEither
    } yield  ()).fold(
      identity,
      _ => original //this should never happen as we're gathering up the errors
    )

  private def decodeCwl(json: Json): ErrorOr[Unit] = {

    val rawJson = io.circe.Printer.noSpaces.pretty(json)

    // we know we are dealing with an individual CWL, so we can take some liberties like assuming "class" is defined.
    ((json \\ "class").head.asString match {
      case Some("Workflow") => decodeWithErrorStrings[Workflow](rawJson)
      case Some("CommandLineTool") => decodeWithErrorStrings[CommandLineTool](rawJson)
      case Some("ExpressionTool") => decodeWithErrorStrings[ExpressionTool](rawJson)
      case _ => "Class field was declared incorrectly!".invalidNel
    }).map(_ => ())
  }
}
