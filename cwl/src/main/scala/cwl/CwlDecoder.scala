package cwl

import ammonite.ops.ImplicitWd._
import ammonite.ops._
import better.files.File.newTemporaryFile
import better.files.{File => BFile}
import cats.Applicative
import cats.data.EitherT._
import cats.data.{EitherT, NonEmptyList, ValidatedNel}
import cats.effect.IO
import cats.syntax.either._
import common.legacy.TwoElevenSupport._
import common.validation.ErrorOr._

import scala.util.Try

object CwlDecoder {

  object Parse {
    def error[A](error: String, tail: String*): Parse[A] = EitherT.leftT {
      NonEmptyList.of(error, tail: _*)
    }
  }
  
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

  def parseJson(json: String): Parse[CwlFile] = fromEither[IO](CwlCodecs.decodeCwl(json))

  /**
   * Notice it gives you one instance of Cwl.  This has transformed all embedded files into scala object state
   */
  def decodeAllCwl(fileName: BFile, root: Option[String] = None): Parse[Cwl] =
    for {
      jsonString <- preprocess(fileName)
      unmodifiedCwl <- parseJson(jsonString)
      cwlWithEmbeddedCwl <- unmodifiedCwl.fold(FlattenCwlFile).apply((fileName.toString, root))
    } yield cwlWithEmbeddedCwl

  def decodeTopLevelCwl(fileName: BFile, rootName: Option[String]): Parse[Cwl] =
    for {
      jsonString <- preprocess(fileName)
      unmodifiedCwl <- parseJson(jsonString)
      rootCwl <- EitherT.fromEither(unmodifiedCwl.fold(FlattenCwlFile.CwlFileRoot).apply(rootName)): Parse[Cwl]
    } yield rootCwl

  def decodeTopLevelCwl(cwl: String, rootName: Option[String] = None): Parse[Cwl] =
    for {
     file <- fromEither[IO](newTemporaryFile().write(cwl).asRight)
     out <- decodeTopLevelCwl(file, rootName)
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

