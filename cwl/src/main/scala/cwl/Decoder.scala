package cwl

import cats.data.{EitherT, NonEmptyList, ValidatedNel}
import ammonite.ops._
import ammonite.ops.ImplicitWd._
import cats.effect.IO
import cats.syntax.either._
import cats.Applicative
import better.files.{File => BFile}
import common.validation.ErrorOr._
import common.legacy.TwoElevenSupport._
import io.circe.{DecodingFailure, ParsingFailure}
import EitherT._
import better.files.File.newTemporaryFile

import scala.util.Try

object CwlDecoder {

  type Parse[A] = EitherT[IO, NonEmptyList[String], A]

  type ParseValidated[A] = IO[ValidatedNel[String, A]]

  implicit val composedApplicative = Applicative[IO] compose Applicative[ErrorOr]


  private def preprocess(path: BFile): Parse[String] = {
    def resultToEither(cr: CommandResult) =
      cr.exitCode match {
        case 0 => Right(cr.out.string)
        case error => Left(NonEmptyList.one(s"running CwlTool on file $path resulted in exit code $error and stderr ${cr.err.string}"))
      }

    val cwlToolResult =
      Try(%%("cwltool", "--quiet", "--print-pre", path.toString)).
        tacticalToEither.
        leftMap(t => NonEmptyList.one(s"running cwltool on file ${path.toString} failed with ${t.getMessage}"))

    fromEither[IO](cwlToolResult flatMap resultToEither)
  }

  def parseJson(json: String): Parse[Cwl] =
    fromEither[IO] {
      CwlCodecs.decodeCwl(json).
        leftMap{
          case df@DecodingFailure(message, ops) => NonEmptyList.of(message, ops.mkString("\n"), df.getStackTrace.mkString("\n"))
          case ParsingFailure(message, underlying) => NonEmptyList.of(message, underlying.getMessage, underlying.getStackTrace.mkString("\n"))
        }
    }

  /**
   * Notice it gives you one instance of Cwl.  This has transformed all embedded files into scala object state
   */
  def decodeAllCwl(fileName: BFile): Parse[Cwl] =
    for {
      jsonString <- preprocess(fileName)
      unmodifiedCwl <- parseJson(jsonString)
      cwlWithEmbeddedCwl <- unmodifiedCwl.fold(AddEmbeddedCwl)
    } yield cwlWithEmbeddedCwl

  def decodeTopLevelCwl(fileName: BFile): Parse[Cwl] =
    for {
      jsonString <- preprocess(fileName)
      unmodifiedCwl <- parseJson(jsonString)
    } yield unmodifiedCwl

  def decodeTopLevelCwl(cwl: String): Parse[Cwl] =
    for {
     file <- fromEither[IO](newTemporaryFile().write(cwl).asRight)
     out <- decodeTopLevelCwl(file)
    } yield out

  //This is used when traversing over Cwl and replacing links w/ embedded data
  private[cwl] def decodeCwlAsValidated(fileName: String): ParseValidated[(String, Cwl)] = {
    //The SALAD preprocess step puts "file://" as a prefix to all filenames.  Better files doesn't like this.
    val bFileName = fileName.drop(5)

    decodeAllCwl(BFile(bFileName)).
      map(fileName.toString -> _).
      value.
      map(_.toValidated)
  }
}

