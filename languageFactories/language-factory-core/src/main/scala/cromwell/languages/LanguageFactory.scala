package cromwell.languages

import common.Checked
import common.validation.Parse.Parse
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.languages.util.ImportResolver.ImportResolver
import wom.core._
import wom.executable.WomBundle
import wom.expression.IoFunctionSet

trait LanguageFactory {

  // Passed in by the constructor:
  def config: Map[String, Any]

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
}

