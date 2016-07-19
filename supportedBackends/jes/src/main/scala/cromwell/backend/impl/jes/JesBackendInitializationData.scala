package cromwell.backend.impl.jes

import com.google.api.services.genomics.Genomics
import cromwell.backend.BackendInitializationData

case class JesBackendInitializationData(workflowPaths: JesWorkflowPaths, genomics: Genomics)
  extends BackendInitializationData
