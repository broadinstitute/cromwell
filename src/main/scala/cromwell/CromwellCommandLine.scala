package cromwell

import java.nio.file.{Files, Path, Paths}

import better.files._
import cromwell.core.{ErrorOr, WorkflowSourceFiles}
import cromwell.util.FileUtil._
import lenthall.exception.MessageAggregation

import scalaz._
import Scalaz._
import scala.util.{Failure, Success, Try}

sealed abstract class CromwellCommandLine
case object UsageAndExit extends CromwellCommandLine
case object RunServer extends CromwellCommandLine
final case class RunSingle(wdlPath: Path,
                           sourceFiles: WorkflowSourceFiles,
                           inputsPath: Option[Path],
                           optionsPath: Option[Path],
                           metadataPath: Option[Path]) extends CromwellCommandLine

object CromwellCommandLine {
  def apply(args: Seq[String]): CromwellCommandLine = {
    args.headOption match {
      case Some("server") if args.size == 1 => RunServer
      case Some("run") if args.size >= 2 && args.size <= 5 => RunSingle(args.tail)
      case _ => UsageAndExit
    }
  }
}

object RunSingle {
  def apply(args: Seq[String]): RunSingle = {
    val wdlPath = Paths.get(args.head)
    val inputsPath = argPath(args, 1, Option(".inputs"), checkDefaultExists = false)
    val optionsPath = argPath(args, 2, Option(".options"), checkDefaultExists = true)
    val metadataPath = argPath(args, 3, None)

    val wdl = readContent("WDL file", wdlPath)
    val inputsJson = readJson("Inputs", inputsPath)
    val optionsJson = readJson("Workflow Options", optionsPath)

    (wdl |@| inputsJson |@| optionsJson) {
      (w, i, o) =>
        val workflowSourceFiles = WorkflowSourceFiles(w, i, o)
        metadataPath foreach { m => if (!metadataPathIsWriteable(m)) throw new RuntimeException(s"Unable to write to metadata directory: $m")}
        RunSingle(wdlPath, workflowSourceFiles, inputsPath, optionsPath, metadataPath)
    } match {
      case scalaz.Success(r) => r
      case scalaz.Failure(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "ERROR: Unable to run Cromwell:"
        override def errorMessages: Traversable[String] = nel.list
      }
    }
  }

  /** Read the path to a string. */
  private def readContent(inputDescription: String, path: Path): ErrorOr[String] = {
    if (!Files.exists(path)) {
      s"$inputDescription does not exist: $path".failureNel
    } else if (!Files.isReadable(path)) {
      s"$inputDescription is not readable: $path".failureNel
    } else path.contentAsString.successNel
  }

  /** Read the path to a string, unless the path is None, in which case returns "{}". */
  private def readJson(inputDescription: String, pathOption: Option[Path]): ErrorOr[String] = {
    pathOption match {
      case Some(path) => readContent(inputDescription, path)
      case None => "{}".successNel
    }
  }

  private def metadataPathIsWriteable(metadataPath: Path): Boolean = {
    Try(metadataPath.createIfNotExists().append("")) match {
      case Success(_) => true
      case Failure(_) => false
    }
  }

  /**
    * Retrieve the arg at index as path, or return some default. Args specified as "-" will be returned as None.
    *
    * @param args The run command arguments, with the wdl path at arg.head.
    * @param index The index of the path we're looking for.
    * @param defaultExt The default extension to use if the argument was not specified at all.
    * @param checkDefaultExists If true, verify that our computed default file exists before using it.
    * @return The argument as a Path resolved as a sibling to the wdl path.
    */
  private def argPath(args: Seq[String], index: Int, defaultExt: Option[String],
                      checkDefaultExists: Boolean = true): Option[Path] = {

    // To return a default, swap the extension, and then maybe check if the file exists.
    def defaultPath = defaultExt
      .map(ext => swapExt(args.head, ".wdl", ext))
      .filter(path => !checkDefaultExists || Files.exists(Paths.get(path)))

    // Return the path for the arg index, or the default, but remove "-" paths.
    for {
      path <- args.lift(index) orElse defaultPath filterNot (_ == "-")
    } yield Paths.get(path)
  }
}
