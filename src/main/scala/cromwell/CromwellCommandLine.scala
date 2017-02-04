package cromwell

import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.{WorkflowSourceFilesCollection, WorkflowSourceFilesWithDependenciesZip, WorkflowSourceFilesWithoutImports}
import lenthall.exception.MessageAggregation
import lenthall.validation.ErrorOr._
import org.slf4j.LoggerFactory

import scala.util.{Failure, Success, Try}

sealed abstract class CromwellCommandLine
case object UsageAndExit extends CromwellCommandLine
case object RunServer extends CromwellCommandLine
case object VersionAndExit extends CromwellCommandLine

object CromwellCommandLine {
  def apply(args: Seq[String]): CromwellCommandLine = {
    args.headOption match {
      case Some("server") if args.size == 1 => RunServer
      case Some("run") if args.size >= 2 && args.size <= 6 => RunSingle(args.tail)
      case Some("-version") if args.size == 1 => VersionAndExit
      case _ => UsageAndExit
    }
  }
}

// We cannot initialize the logging until after we parse the command line in Main.scala. So we have to bundle up and pass back this information, just for logging.
case class SingleRunPathParameters(wdlPath: Path, inputsPath: Option[Path], optionsPath: Option[Path], metadataPath: Option[Path], importPath: Option[Path], labelsPath: Option[Path]) {
  def logMe(log: org.slf4j.Logger) = {
    log.info(s"  WDL file: $wdlPath")
    inputsPath foreach { i => log.info(s"  Inputs: $i") }
    optionsPath foreach { o => log.info(s"  Workflow Options: $o") }
    metadataPath foreach { m => log.info(s"  Workflow Metadata Output: $m") }
    importPath foreach { i => log.info(s"  WDL import bundle: $i") }
    labelsPath foreach { o => log.info(s"  Custom labels: $o") }
  }
}

final case class RunSingle(sourceFiles: WorkflowSourceFilesCollection, paths: SingleRunPathParameters) extends CromwellCommandLine

object RunSingle {

  lazy val Log = LoggerFactory.getLogger("cromwell")

  def apply(args: Seq[String]): RunSingle = {
    val pathParameters = SingleRunPathParameters(
      wdlPath = DefaultPathBuilder.get(args.head).toAbsolutePath,
      inputsPath = argPath(args, 1, Option("inputs"), checkDefaultExists = false),
      optionsPath = argPath(args, 2, Option("options")),
      metadataPath = argPath(args, 3, None),
      importPath = argPath(args, 4, None),
      labelsPath = argPath(args, 5, None)
    )

    val wdl = readContent("WDL file", pathParameters.wdlPath)
    val inputsJson = readJson("Inputs", pathParameters.inputsPath)
    val optionsJson = readJson("Workflow Options", pathParameters.optionsPath)
    val labelsJson = readJson("Labels", pathParameters.labelsPath)

    val sourceFileCollection = pathParameters.importPath match {
      case Some(p) => (wdl |@| inputsJson |@| optionsJson |@| labelsJson) map { (w, i, o, l) =>
        WorkflowSourceFilesWithDependenciesZip.apply(w, i, o, l, p.loadBytes) }
      case None => (wdl |@| inputsJson |@| optionsJson |@| labelsJson) map WorkflowSourceFilesWithoutImports.apply
    }

    val runSingle: ErrorOr[RunSingle] = for {
      sources <- sourceFileCollection
      _ <- writeableMetadataPath(pathParameters.metadataPath)
    } yield RunSingle(sources, pathParameters)

    runSingle match {
      case Valid(r) => r
      case Invalid(nel) => throw new RuntimeException with MessageAggregation {
        override def exceptionContext: String = "ERROR: Unable to run Cromwell:"
        override def errorMessages: Traversable[String] = nel.toList
      }
    }
  }

  private def writeableMetadataPath(path: Option[Path]): ErrorOr[Unit] = {
    path match {
      case Some(p) if !metadataPathIsWriteable(p) => s"Unable to write to metadata directory: $p".invalidNel
      case _ => ().validNel
    }
  }

  /** Read the path to a string. */
  private def readContent(inputDescription: String, path: Path): ErrorOr[String] = {
    if (!path.exists) {
      s"$inputDescription does not exist: $path".invalidNel
    } else if (!path.isReadable) {
      s"$inputDescription is not readable: $path".invalidNel
    } else path.contentAsString.validNel
  }

  /** Read the path to a string, unless the path is None, in which case returns "{}". */
  private def readJson(inputDescription: String, pathOption: Option[Path]): ErrorOr[String] = {
    pathOption match {
      case Some(path) => readContent(inputDescription, path)
      case None => "{}".validNel
    }
  }

  private def metadataPathIsWriteable(metadataPath: Path): Boolean = {
    Try(metadataPath.createIfNotExists(createParents = true).append("")) match {
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
      .map(ext => DefaultPathBuilder.get(args.head).swapExt("wdl", ext))
      .filter(path => !checkDefaultExists || path.exists)
      .map(_.pathAsString)

    // Return the path for the arg index, or the default, but remove "-" paths.
    for {
      path <- args.lift(index) orElse defaultPath filterNot (_ == "-")
    } yield DefaultPathBuilder.get(path).toAbsolutePath
  }
}
