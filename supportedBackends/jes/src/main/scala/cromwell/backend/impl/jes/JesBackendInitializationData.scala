package cromwell.backend.impl.jes

import com.google.api.services.genomics.Genomics
import cromwell.backend.standard.StandardInitializationData

case class JesBackendInitializationData
(
  override val workflowPaths: JesWorkflowPaths,
  jesConfiguration: JesConfiguration,
  genomics: Genomics
) extends StandardInitializationData(workflowPaths, JesRuntimeAttributes.runtimeAttributesBuilder)
