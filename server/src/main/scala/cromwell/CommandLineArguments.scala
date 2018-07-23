package cromwell

import java.net.URL

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import cats.syntax.apply._
import cats.syntax.either._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Parse._
import common.validation.Validation._
import cromwell.CommandLineArguments._
import cromwell.CromwellApp.Command
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.webservice.CromwellApiService
import cwl.preprocessor.CwlPreProcessor
import org.slf4j.Logger

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

object CommandLineArguments {
  val DefaultCromwellHost = new URL("http://localhost:8000")
  case class ValidSubmission(
                              workflowSource: String,
                              workflowUrl: Option[String],
                              workflowRoot: Option[String],
                              worflowInputs: String,
                              workflowOptions: String,
                              workflowLabels: String,
                              dependencies: Option[File])
}

case class CommandLineArguments(command: Option[Command] = None,
                                workflowSource: Option[Path] = None,
                                workflowUrl: Option[String] = None,
                                workflowRoot: Option[String] = None,
                                workflowInputs: Option[Path] = None,
                                workflowOptions: Option[Path] = None,
                                workflowType: Option[String] = None,
                                workflowTypeVersion: Option[String] = None,
                                workflowLabels: Option[Path] = None,
                                imports: Option[Path] = None,
                                metadataOutput: Option[Path] = None,
                                host: URL = CommandLineArguments.DefaultCromwellHost
                               ) {
  private lazy val cwlPreProcessor = new CwlPreProcessor()
  private lazy val isCwl = workflowType.exists(_.equalsIgnoreCase("cwl"))

  /*
    * If the file is a relative local path, resolve it against the path of the input json.
   */
  private def inputFilesMapper(inputJsonPath: Path)(file: String) = {
    DefaultPathBuilder.build(file) match {
      case Success(path) if !path.isAbsolute => inputJsonPath.sibling(file).pathAsString
      case _ => file
    }
  }

  private def preProcessCwlInputFile(path: Path): ErrorOr[String] = {
    cwlPreProcessor.preProcessInputFiles(path.contentAsString, inputFilesMapper(path)).toErrorOr
  }

  def validateSubmission(logger: Logger)
                        (implicit ec: ExecutionContext, materializer: ActorMaterializer, actorSystem: ActorSystem): Future[ErrorOr[ValidSubmission]] = {
    val workflowPath = File(workflowSource.get.pathAsString)

    val urlContentFuture: Future[Option[String]] = CromwellApiService.getContentFromWorkflowUrl(workflowUrl)(ec, materializer, actorSystem)

    urlContentFuture map { urlContent =>

      val workflowAndDependencies: ErrorOr[(String, Option[File], Option[String])] = if (isCwl) {
        logger.info("Pre Processing Workflow...")
        lazy val preProcessedCwl = cwlPreProcessor.preProcessCwlFileToString(workflowPath, None)

        imports match {
          case Some(explicitImports) => readOptionContent("Workflow source", workflowSource).map((_, Option(File(explicitImports.pathAsString)), workflowRoot))
          case None => Try(preProcessedCwl.map((_, None, None)).value.unsafeRunSync())
            .toChecked
            .flatMap(identity)
            .toValidated
        }
      } else readOptionContent("Workflow source", workflowSource).map((_, imports.map(p => File(p.pathAsString)), workflowRoot))

      val inputsJson: ErrorOr[String] = if (isCwl) {
        logger.info("Pre Processing Inputs...")
        workflowInputs.map(preProcessCwlInputFile).getOrElse(readOptionContent("Workflow inputs", workflowInputs))
      } else readOptionContent("Workflow inputs", workflowInputs)

      val optionsJson = readOptionContent("Workflow options", workflowOptions)
      val labelsJson = readOptionContent("Workflow labels", workflowLabels)

      val workflowSourceFinal: ErrorOr[String] = (workflowSource, workflowUrl) match {
        case (Some(path), None) => readContent("Workflow source", path)
        case (None, Some(_)) => "".validNel //urlContent
        case (Some(_), Some(_)) => "Both Workflow source and Workflow url can't be supplied".invalidNel
        case (None, None) => "Workflow source and Workflow url needs to be supplied".invalidNel
      }

      (workflowAndDependencies, inputsJson, optionsJson, labelsJson, workflowSourceFinal) mapN {
        case ((_, z, r), i, o, l, w) =>
          ValidSubmission(w, workflowUrl, r, i, o, l, z)
      }
    }
  }

//  def validateSubmissionCopy(logger: Logger, urlContent: Option[String]): ErrorOr[ValidSubmission] = {
//    val workflowPath = File(workflowSource.get.pathAsString)
//
//    val workflowAndDependencies: ErrorOr[(String, Option[File], Option[String])] = if (isCwl) {
//      logger.info("Pre Processing Workflow...")
//      lazy val preProcessedCwl = cwlPreProcessor.preProcessCwlFileToString(workflowPath, None)
//
//      imports match {
//        case Some(explicitImports) => readOptionContent("Workflow source", workflowSource).map((_, Option(File(explicitImports.pathAsString)), workflowRoot))
//        case None => Try(preProcessedCwl.map((_, None, None)).value.unsafeRunSync())
//          .toChecked
//          .flatMap(identity)
//          .toValidated
//      }
//    } else readOptionContent("Workflow source", workflowSource).map((_, imports.map(p => File(p.pathAsString)), workflowRoot))
//
//    val inputsJson: ErrorOr[String] = if (isCwl) {
//      logger.info("Pre Processing Inputs...")
//      workflowInputs.map(preProcessCwlInputFile).getOrElse(readOptionContent("Workflow inputs", workflowInputs))
//    } else readOptionContent("Workflow inputs", workflowInputs)
//
//    val optionsJson = readOptionContent("Workflow options", workflowOptions)
//    val labelsJson = readOptionContent("Workflow labels", workflowLabels)
//
//    val workflowSourceFinal: ErrorOr[String] = (workflowSource, workflowUrl) match {
//      case (Some(path), None) => readContent("Workflow source", path)
//      case (_, Some(url)) => url.validNel
//      case (Some(_), Some(_)) => "Workflow source and Workflow url can't both be supplied".invalidNel
//      case (None, None) => "Workflow source and Workflow url needs to be supplied".invalidNel
//    }
//
//    //    for {
//    //      _ <- CromwellApiService.downloadContentFromUrl(workflowUrl)
//    //    } yield ()
//
//    (workflowAndDependencies, inputsJson, optionsJson, labelsJson, workflowSourceFinal) mapN {
//      case ((w, z, r), i, o, l, _) =>
//        ValidSubmission(w, workflowUrl, r, i, o, l, z)
//    }
//  }

//  def abc(logger: Logger): ErrorOr[ValidSubmission] = {
//    val urlContentFuture = CromwellApiService.getContentFromWorkflowUrl(workflowUrl)
//    val urlContent = Await.result(urlContentFuture, Duration.Inf)

    for{
      urlContent <- CromwellApiService.getContentFromWorkflowUrl(workflowUrl)
      submission <- Future.successful(validateSubmissionCopy(logger, urlContent))
    } yield submission
//
//    validateSubmissionCopy(logger, urlContent)
//  }

  /** Read the path to a string. */
  private def readContent(inputDescription: String, path: Path): ErrorOr[String] = {
    if (!path.exists) {
      s"$inputDescription does not exist: $path".invalidNel
    } else if (!path.isReadable) {
      s"$inputDescription is not readable: $path".invalidNel
    } else path.contentAsString.validNel
  }

  /** Read the path to a string, unless the path is None, in which case returns "{}". */
  private def readOptionContent(inputDescription: String, pathOption: Option[Path]): ErrorOr[String] = {
    pathOption match {
      case Some(path) => readContent(inputDescription, path)
      case None => "{}".validNel
    }
  }
}
