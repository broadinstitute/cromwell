package womtool.input

import java.nio.file.{Files, Paths}

import com.typesafe.config.ConfigFactory
import common.Checked
import common.validation.Validation._
import cromwell.core.path.Path
import cromwell.languages.LanguageFactory
import cromwell.languages.util.ImportResolver._
import languages.cwl.CwlV1_0LanguageFactory
import languages.wdl.biscayne.WdlBiscayneLanguageFactory
import languages.wdl.draft2.WdlDraft2LanguageFactory
import languages.wdl.draft3.WdlDraft3LanguageFactory
import wom.ResolvedImportRecord
import wom.executable.WomBundle
import wom.expression.NoIoFunctionSet
import wom.graph._

import scala.collection.JavaConverters._
import scala.util.Try

object WomGraphMaker {

  def getBundle(mainFile: Path): Checked[WomBundle] = getBundleAndFactory(mainFile).map(_._1)

  private def getBundleAndFactory(mainFile: Path): Checked[(WomBundle, LanguageFactory)] = {
    lazy val importResolvers: List[ImportResolver] =
      DirectoryResolver.localFilesystemResolvers(Some(mainFile)) :+ HttpResolver(relativeTo = None)

    readFile(mainFile.toAbsolutePath.pathAsString) flatMap { mainFileContents =>
      val languageFactory =
        List(
          new WdlDraft3LanguageFactory(ConfigFactory.empty()),
          new WdlBiscayneLanguageFactory(ConfigFactory.empty()),
          new CwlV1_0LanguageFactory(ConfigFactory.empty()))
          .find(_.looksParsable(mainFileContents))
          .getOrElse(new WdlDraft2LanguageFactory(ConfigFactory.empty()))

      val bundle = languageFactory.getWomBundle(mainFileContents, None, "{}", importResolvers, List(languageFactory))
      // Return the pair with the languageFactory
      bundle map ((_, languageFactory))
    }
  }

  def fromFiles(mainFile: Path, inputs: Option[Path]): Checked[WomGraphWithResolvedImports] = {
    getBundleAndFactory(mainFile) flatMap { case (womBundle, languageFactory) =>
      inputs match {
        case None =>
          for {
            executableCallable <- womBundle.toExecutableCallable
          } yield WomGraphWithResolvedImports(executableCallable.graph, womBundle.resolvedImportRecords)
        case Some(inputsFile) =>
          for {
            inputsContents <- readFile(inputsFile.toAbsolutePath.pathAsString)
            validatedWomNamespace <- languageFactory.createExecutable(womBundle, inputsContents, NoIoFunctionSet)
          } yield WomGraphWithResolvedImports(validatedWomNamespace.executable.graph, womBundle.resolvedImportRecords)
      }
    }
  }

  private def readFile(filePath: String): Checked[String] = Try(Files.readAllLines(Paths.get(filePath)).asScala.mkString(System.lineSeparator())).toChecked

}


case class WomGraphWithResolvedImports(graph: Graph, resolvedImportRecords: Set[ResolvedImportRecord])
