package cromwell.backend.impl

import com.google.api.services.genomics.Genomics
import cromwell.backend.BackendInitializationData

package object jes {
  case class JesBackendInitializationData(workflowPaths: JesWorkflowPaths, genomics: Genomics) extends BackendInitializationData
}
