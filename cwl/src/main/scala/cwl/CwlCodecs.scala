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
import io.circe.Json._
import io.circe.ParsingFailure._
import io.circe.DecodingFailure._

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

  /**
    * Attempt to parse as an array of CWL or a single one.
    *
    * If it fails, attempt to parse each CWL individually and return all of those failures.
    */
  def decodeCwl(in: String): Checked[CwlFile] =
    decodeWithErrorStrings[CwlFile](in).leftMap(_ => decodePieces(in)).toEither

  private def decodeWithErrorStrings[A](in: String)(implicit d: Decoder[A]): ErrorOr[A] =
    decodeAccumulating[A](in).leftMap(_.map(_.show))

  private def decodeWithErrorStringsJson[A](in: Json)(implicit d: Decoder[A]): ErrorOr[A] =
    in.as[A].leftMap(_.show).leftMap(NonEmptyList.one).toValidated

  /**
    * Drop down to a lower level and try to parse each CWL individually.
    */
  private def decodePieces(originalJson: String):NonEmptyList[String]  = {

    //break down a CWL into pieces, and try to decode each one.  Gather up all the errors
    val errorsPerPiece: Either[NonEmptyList[String], Unit] =
      for {
        raw <- parse(originalJson).leftMap(_.show).leftMap(NonEmptyList.one)
        _ <-
          if (raw.isArray)
            raw.asArray.toVector.flatten.traverse(decodeCwl).toEither
          else
            decodeCwl(raw).toEither
      } yield ()

    //we need a way to get the errors and ignore the "success" case. Validated doesn't have a "getOrElse" for errors so we're stuck with fold.
    errorsPerPiece.fold(
      identity,
      _ => throw new RuntimeException("There is a bug in Cromwell's parsing logic.  The code meant to improve error messaging has encountered an unknown state")
    )
  }

  private def findClass(json: Json): Option[String] =
    for {
      obj <- json.asObject
      map = obj.toMap
      classObj <- map.get("class")
      classString <- classObj.asString
    } yield classString

  private def decodeCwl(json: Json): ErrorOr[Unit] = {
      (findClass(json) match {
      case Some("Workflow") => decodeWithErrorStringsJson[Workflow](json)
      case Some("CommandLineTool") => decodeWithErrorStringsJson[CommandLineTool](json)
      case Some("ExpressionTool") => decodeWithErrorStringsJson[ExpressionTool](json)
      case Some(other) => s"Class field was declared incorrectly: $other is not one of Workflow, CommandLineTool, or ExpressionTool! as seen in ${json.show}".invalidNel
      case None => s"Class field was omitted in ${json.show}".invalidNel
    }).map(_ => ())
  }
}
