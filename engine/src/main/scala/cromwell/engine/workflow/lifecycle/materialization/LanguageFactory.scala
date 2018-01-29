package cromwell.engine.workflow.lifecycle.materialization

import akka.actor.ActorRef
import common.validation.Parse.Parse
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesCollection}
import wom.executable.ValidatedWomNamespace
import wom.expression.IoFunctionSet

trait LanguageFactory {
  def validateNamespace(source: WorkflowSourceFilesCollection,
                           workflowOptions: WorkflowOptions,
                           ioFunctions: IoFunctionSet,
                           importLocalFilesystem: Boolean,
                           serviceRegistryActor: ActorRef,
                           workflowIdForLogging: WorkflowId): Parse[ValidatedWomNamespace]
}
