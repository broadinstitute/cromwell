package cromwell.backend.impl.local

import cromwell.backend.io.{WorkflowPaths, WorkflowPathsBackendInitializationData}

case class LocalBackendInitializationData(workflowPaths: WorkflowPaths) extends WorkflowPathsBackendInitializationData
