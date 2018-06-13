package womtool.input

import java.nio.file.{Files, Paths}

import common.Checked
import common.validation.Validation._
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.languages.LanguageFactory
import cromwell.languages.util.ImportResolver._
import languages.cwl.CwlV1_0LanguageFactory
import languages.wdl.draft2.WdlDraft2LanguageFactory
import languages.wdl.draft3.WdlDraft3LanguageFactory
import wom.executable.WomBundle
import wom.expression.NoIoFunctionSet
import wom.graph._

import scala.collection.JavaConverters._
import scala.util.Try

object WomGraphMaker {

  def getBundle(mainFile: Path): Checked[WomBundle] = getBundleAndFactory(mainFile).map(_._1)

  private def getBundleAndFactory(mainFile: Path): Checked[(WomBundle, LanguageFactory)] = {
    // Resolves for:
    // - Where we run from
    // - Where the file is
    lazy val importResolvers = List(
      directoryResolver(DefaultPathBuilder.build(Paths.get("."))),
      directoryResolver(DefaultPathBuilder.build(Paths.get(mainFile.toAbsolutePath.toFile.getParent)))
    )

    readFile(mainFile.toAbsolutePath.pathAsString) flatMap { mainFileContents =>
      val languageFactory = if (mainFile.name.toLowerCase().endsWith("wdl")) {
        if (mainFileContents.startsWith("version 1.0") || mainFileContents.startsWith("version draft-3")) {
          new WdlDraft3LanguageFactory(Map.empty)
        } else {
          new WdlDraft2LanguageFactory(Map.empty)
        }
      } else new CwlV1_0LanguageFactory(Map.empty)

      val bundle = languageFactory.getWomBundle(mainFileContents, "{}", importResolvers, List(languageFactory))
      // Return the pair with the languageFactory
      bundle map ((_, languageFactory))
    }
  }

  def fromFiles(mainFile: Path, inputs: Option[Path]): Checked[Graph] = {
    getBundleAndFactory(mainFile) flatMap { case (womBundle, languageFactory) =>
      inputs match {
        case None =>
          for {
            executableCallable <- womBundle.toExecutableCallable
          } yield executableCallable.graph
        case Some(inputsFile) =>
          for {
            inputsContents <- readFile(inputsFile.toAbsolutePath.pathAsString)
            validatedWomNamespace <- languageFactory.createExecutable(womBundle, inputsContents, NoIoFunctionSet)
          } yield validatedWomNamespace.executable.graph
      }
    }
  }

  private def readFile(filePath: String): Checked[String] = Try(Files.readAllLines(Paths.get(filePath)).asScala.mkString(System.lineSeparator())).toChecked

}
