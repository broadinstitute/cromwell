package cromwell.languages

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
  def config: Map[String, Any]

  def enabledCheck: Checked[Unit] = if (standardConfig.enabled) {
    ().validNelCheck
  } else {
    s"The language factory for $languageName ($languageVersionName) is not currently enabled in this Cromwell".invalidNelCheck
  }

  // Override if you want to accumulate extra options instead of throwing an exception:
  lazy val standardConfig: StandardLanguageFactoryConfig = StandardLanguageFactoryConfig.parse(config, allowExtras = false)

  def getWomBundle(workflowSource: WorkflowSource,
                   workflowOptionsJson: WorkflowOptionsJson,
                   importResolvers: List[ImportResolver],
                   languageFactories: List[LanguageFactory]): Checked[WomBundle]

  def createExecutable(womBundle: WomBundle,
                       inputs: WorkflowJson,
                       ioFunctions: IoFunctionSet): Checked[ValidatedWomNamespace]

  def validateNamespace(source: WorkflowSourceFilesCollection,
                        workflowOptions: WorkflowOptions,
                        importLocalFilesystem: Boolean,
                        workflowIdForLogging: WorkflowId,
                        ioFunctions: IoFunctionSet): Parse[ValidatedWomNamespace]

  /**
    * In case no version is specified: does this language factory feel like it might be suitable for this file?
    * @param content The workflow description
    */
  def looksParsable(content: String): Boolean
}
