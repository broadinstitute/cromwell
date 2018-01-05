package cromwell.backend.impl.bcs

import cromwell.backend.io.WorkflowPaths
import cromwell.backend.standard.{StandardInitializationData, StandardValidatedRuntimeAttributesBuilder}
import cromwell.core.path.PathBuilder

final case class BcsBackendInitializationData
(
  override val workflowPaths: WorkflowPaths,
  override val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder,
  bcsConfiguration: BcsConfiguration,
  pathBuilders: List[PathBuilder]
) extends StandardInitializationData(workflowPaths, runtimeAttributesBuilder, classOf[BcsExpressionFunctions])

