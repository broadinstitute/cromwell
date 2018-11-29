package cwl

import better.files.{File => BFile}
import cats.data.EitherT._
import cats.effect.IO
import cats.syntax.either._
import cats.{Applicative, Monad}
import common.validation.ErrorOr._
import common.validation.IOChecked._
import common.validation.Validation._
import cwl.preprocessor.{CwlFileReference, CwlPreProcessor, CwlReference}
import io.circe.Json

import scala.util.Try

object CwlDecoder {

  implicit val composedApplicative = Applicative[IO] compose Applicative[ErrorOr]

  def saladCwlFile(reference: CwlReference): IOChecked[String] = {
    val cwlToolResult =
      Try(CwltoolRunner.instance.salad(reference))
        .toCheckedWithContext(s"run cwltool on file ${reference.pathAsString}", throwableToStringFunction = t => t.toString)

    fromEither[IO](cwlToolResult)
  }

  private lazy val cwlPreProcessor = new CwlPreProcessor()

  // TODO: WOM: During conformance testing the saladed-CWLs are referring to files in the temp directory.
  // Thus we can't delete the temp directory until after the workflow is complete, like the workflow logs.
  // All callers to this method should be fixed around the same time.
  // https://github.com/broadinstitute/cromwell/issues/3186
  def todoDeleteCwlFileParentDirectory(cwlFile: BFile): IOChecked[Unit] = {
    goIOChecked {
      //cwlFile.parent.delete(swallowIOExceptions = true)
    }
  }

  def parseJson(json: Json, from: String): IOChecked[Cwl] = fromEither[IO](CwlCodecs.decodeCwl(json).contextualizeErrors(s"parse '$from'"))

  /**
    * Notice it gives you one instance of Cwl.  This has transformed all embedded files into scala object state
    */
  def decodeCwlReference(reference: CwlReference)(implicit processor: CwlPreProcessor = cwlPreProcessor): IOChecked[Cwl] = {
    def makeStandaloneWorkflow(): IOChecked[Json] = processor.preProcessCwl(reference)

    for {
      standaloneWorkflow <- makeStandaloneWorkflow()
      parsedCwl <- parseJson(standaloneWorkflow, reference.pathAsString)
    } yield parsedCwl
  }

  def decodeCwlFile(file: BFile, workflowRoot: Option[String] = None) = {
    decodeCwlReference(CwlFileReference(file, workflowRoot))
  }

  def decodeCwlString(cwl: String,
                      zipOption: Option[BFile] = None,
                      rootName: Option[String] = None,
                      cwlFilename: String = "cwl_temp_file"): IOChecked[Cwl] = {
    for {
      parentDir <- goIOChecked(BFile.newTemporaryDirectory("cwl_temp_dir_")) // has a random long appended like `cwl_temp_dir_100000000000`
      file <- fromEither[IO](parentDir./(cwlFilename + ".cwl").write(cwl).asRight) // serves as the basis for the output directory name; must remain stable across restarts
      _ <- zipOption match {
        case Some(zip) => goIOChecked(zip.unzipTo(parentDir))
        case None => Monad[IOChecked].unit
      }
      out <- decodeCwlFile(file, rootName)
      _ <- todoDeleteCwlFileParentDirectory(file)
    } yield out
  }

  //This is used when traversing over Cwl and replacing links w/ embedded data
  private[cwl] def decodeCwlAsValidated(fileName: String): IOCheckedValidated[(String, Cwl)] = {
    //The SALAD preprocess step puts "file://" as a prefix to all filenames.  Better files doesn't like this.
    val bFileName = fileName.stripPrefix("file://")

    decodeCwlReference(CwlFileReference(BFile(bFileName), None)).
      map(fileName.toString -> _).
      value.
      map(_.toValidated)
  }
}
