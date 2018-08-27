package cwl

import better.files.{File => BFile}
import cats.data.EitherT._
import cats.data.NonEmptyList
import cats.effect.IO
import cats.syntax.either._
import cats.{Applicative, Monad}
import common.validation.ErrorOr._
import common.validation.Parse._
import cwl.preprocessor.CwlPreProcessor
import io.circe.Json

import scala.util.Try

object CwlDecoder {

  implicit val composedApplicative = Applicative[IO] compose Applicative[ErrorOr]

  def saladCwlFile(path: BFile): Parse[String] = {
    val cwlToolResult =
      Try(CwltoolRunner.instance.salad(path))
        .toEither
        .leftMap(t => NonEmptyList.one(s"running cwltool on file ${path.toString} failed with ${t.getMessage}"))

    fromEither[IO](cwlToolResult)
  }

  private lazy val cwlPreProcessor = new CwlPreProcessor()

  // TODO: WOM: During conformance testing the saladed-CWLs are referring to files in the temp directory.
  // Thus we can't delete the temp directory until after the workflow is complete, like the workflow logs.
  // All callers to this method should be fixed around the same time.
  // https://github.com/broadinstitute/cromwell/issues/3186
  def todoDeleteCwlFileParentDirectory(cwlFile: BFile): Parse[Unit] = {
    goParse {
      //cwlFile.parent.delete(swallowIOExceptions = true)
    }
  }

  def parseJson(json: Json, file: BFile): Parse[Cwl] = fromEither[IO](CwlCodecs.decodeCwl(json).leftMap(_.prepend(s"error when parsing file $file")))

  /**
    * Notice it gives you one instance of Cwl.  This has transformed all embedded files into scala object state
    */
  def decodeCwlFile(fileName: BFile,
                    workflowRoot: Option[String] = None)(implicit processor: CwlPreProcessor = cwlPreProcessor): Parse[Cwl] =
    for {
      standaloneWorkflow <- processor.preProcessCwlFile(fileName, workflowRoot)
      parsedCwl <- parseJson(standaloneWorkflow, fileName)
    } yield parsedCwl

  def decodeCwlString(cwl: String, zipOption: Option[BFile] = None, rootName: Option[String] = None): Parse[Cwl] = {
    for {
      parentDir <- goParse(BFile.newTemporaryDirectory("cwl.temp."))
      file <- fromEither[IO](BFile.newTemporaryFile("temp.", ".cwl", Option(parentDir)).write(cwl).asRight)
      _ <- zipOption match {
        case Some(zip) => goParse(zip.unzipTo(parentDir))
        case None => Monad[Parse].unit
      }
      out <- decodeCwlFile(file, rootName)
      _ <- todoDeleteCwlFileParentDirectory(file)
    } yield out
  }

  //This is used when traversing over Cwl and replacing links w/ embedded data
  private[cwl] def decodeCwlAsValidated(fileName: String): ParseValidated[(String, Cwl)] = {
    //The SALAD preprocess step puts "file://" as a prefix to all filenames.  Better files doesn't like this.
    val bFileName = fileName.stripPrefix("file://")

    decodeCwlFile(BFile(bFileName)).
      map(fileName.toString -> _).
      value.
      map(_.toValidated)
  }
}

