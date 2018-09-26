package cromwell.languages

import com.typesafe.config.Config
import common.Checked
import common.validation.Parse.Parse
import common.validation.Checked._
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.languages.util.ImportResolver.ImportResolver
import wom.core._
import wom.executable.WomBundle
import wom.expression.IoFunctionSet

trait LanguageFactory {

  def languageName: String
  def languageVersionName: String

  // Passed in by the constructor:
  def config: Config

  import net.ceedubs.ficus.Ficus._

  lazy val enabled = !config.as[Option[Boolean]]("enabled").contains(false)
  lazy val enabledCheck: Checked[Unit] = if (enabled) ().validNelCheck else
    s"The language factory for $languageName ($languageVersionName) is not currently enabled in this Cromwell".invalidNelCheck


  lazy val strictValidation: Boolean = !config.as[Option[Boolean]]("strict-validation").contains(false)

  def getWomBundle(workflowSource: WorkflowSource,
                   workflowOptionsJson: WorkflowOptionsJson,
                   importResolvers: List[ImportResolver],
                   languageFactories: List[LanguageFactory]): Checked[WomBundle]

  def createExecutable(womBundle: WomBundle,
                       inputs: WorkflowJson,
                       ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace]

  def validateNamespace(source: WorkflowSourceFilesCollection,
                        workflowSource: WorkflowSource,
                        workflowOptions: WorkflowOptions,
                        importLocalFilesystem: Boolean,
                        workflowIdForLogging: WorkflowId,
                        ioFunctions: IoFunctionSet,
                        importResolvers: List[ImportResolver]): Parse[ValidatedWomNamespace]

  /**
    * In case no version is specified: does this language factory feel like it might be suitable for this file?
    * @param content The workflow description
    */
  def looksParsable(content: String): Boolean
}
