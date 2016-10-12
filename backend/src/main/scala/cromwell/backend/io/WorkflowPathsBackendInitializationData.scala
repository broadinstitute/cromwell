package cromwell.backend.io

import cromwell.backend.BackendInitializationData
import cromwell.core.path.PathBuilder

/**
  * Extension of backend initialization data that also provides a `WorkflowPaths`, and by proxy its `List[FileSystem]`.
  *
  * This class consolidates code used by the legacy Local backend and the Shared File System (SFS) backend, that need to
  * localize data from a collection of different file systems (local and GCS for now).
  *
  * NOTE: The JES backend uses a fork of the `WorkflowPaths`. As `JesWorkflowPaths` is not a subclass, the
  * initialization data for JES does not use this trait nor the companion object. It is not clear if JES needs the
  * functionality, as it cannot localize from generic file systems, only GCS with different authentication modes.
  *
  * Meanwhile the local and SFS backends need access to the functionality in workflow paths. Their initialization data
  * will therefore be a subset of this trait.
  *
  * The initialization data will be generated in the intialization actor as class extending this trait.
  *
  * From there, the initialization actor will pass the data back to the backend factory, wrapped in an Option, with its
  * type downcast to an `Option[BackendInitializationData]`.
  *
  * The factory then sends a copy of the intialization data down to each execution actor created for every `wdl4s.Call`.
  *
  * Each instance of the execution actor will unpack the initialization data, obtaining the workflow paths using
  * `WorkflowPathsBackendInitializationData.workflowPaths(Option[BackendInitializationData]): WorkflowPaths`
  *
  * However, the execution instances are most interested in actually getting to the file systems embedded within the
  * workflow paths. That functionality is provided via `WorkflowPathsBackendInitializationData.workflowPaths`
  */
trait WorkflowPathsBackendInitializationData extends BackendInitializationData {
  def workflowPaths: WorkflowPaths
}

object WorkflowPathsBackendInitializationData {
  def workflowPaths(initializationDataOption: Option[BackendInitializationData]): WorkflowPaths = {
    BackendInitializationData.as[WorkflowPathsBackendInitializationData](initializationDataOption).workflowPaths
  }

  def pathBuilders(initializationDataOption: Option[BackendInitializationData]): List[PathBuilder] = {
    workflowPaths(initializationDataOption).pathBuilders
  }
}
