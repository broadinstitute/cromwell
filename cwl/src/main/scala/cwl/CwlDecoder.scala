package cwl

import ammonite.ops.ImplicitWd._
import ammonite.ops._
import better.files.{File => BFile}
import cats.data.EitherT._
import cats.data.{EitherT, NonEmptyList, ValidatedNel}
import cats.effect.IO
import cats.syntax.either._
import cats.{Applicative, Monad}
import common.legacy.TwoElevenSupport._
import common.validation.ErrorOr._
import common.validation.Validation._

import scala.util.Try

object CwlDecoder {

  object Parse {
    def error[A](error: String, tail: String*): Parse[A] = EitherT.leftT {
      NonEmptyList.of(error, tail: _*)
    }
  }
  
  type Parse[A] = EitherT[IO, NonEmptyList[String], A]

  type ParseValidated[A] = IO[ValidatedNel[String, A]]

  // If anyone has a magic import that does exactly what these helpers do, please replace, thx!

  def errorOrParse[A](f: => ErrorOr[A]): Parse[A] = fromEither[IO](f.toEither)

  def tryParse[A](f: => Try[A]): Parse[A] = errorOrParse(f.toErrorOr)

  def goParse[A](f: => A): Parse[A] = tryParse(Try(f))

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

  // TODO: WOM: During conformance testing the saladed-CWLs are referring to files in the temp directory.
  // Thus we can't delete the temp directory until after the workflow is complete, like the workflow logs.
  // All callers to this method should be fixed around the same time.
  // https://github.com/broadinstitute/cromwell/issues/3186
  def todoDeleteCwlFileParentDirectory(cwlFile: BFile): Parse[Unit] = {
    goParse {
      //cwlFile.parent.delete(swallowIOExceptions = true)
    }
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

  def decodeTopLevelCwl(cwl: String, zipOption: Option[BFile] = None, rootName: Option[String] = None): Parse[Cwl] = {
    for {
      parentDir <- goParse(BFile.newTemporaryDirectory("cwl.temp."))
      file <- fromEither[IO](BFile.newTemporaryFile("temp.", ".cwl", Option(parentDir)).write(cwl).asRight)
      _ <- zipOption match {
        case Some(zip) => goParse(zip.unzipTo(parentDir))
        case None => Monad[Parse].unit
      }
      out <- decodeTopLevelCwl(file, rootName)
      _ <- todoDeleteCwlFileParentDirectory(file)
    } yield out
  }

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

