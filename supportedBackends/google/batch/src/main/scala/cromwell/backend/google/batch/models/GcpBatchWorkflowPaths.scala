package cromwell.backend.google.batch.models

import com.google.auth.Credentials
import com.typesafe.config.Config
import cromwell.backend.google.batch.runnable.WorkflowOptionKeys
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.WorkflowOptions
import cromwell.core.path.Path
import cromwell.core.path.PathFactory.PathBuilders

import scala.concurrent.ExecutionContext

object GcpBatchWorkflowPaths {
  val GcsRootOptionKey = "jes_gcs_root"
}
case class GcpBatchWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor,
                                 gcsCredentials: Credentials,
                                 batchCredentials: Credentials,
                                 gcpBatchConfiguration: GcpBatchConfiguration,
                                 override val pathBuilders: PathBuilders,
                                 standardStreamNameToFileNameMetadataMapper: (GcpBatchJobPaths, String) => String
)(implicit ec: ExecutionContext)
    extends WorkflowPaths {

  override lazy val executionRootString: String =
    workflowDescriptor.workflowOptions.getOrElse(GcpBatchWorkflowPaths.GcsRootOptionKey, gcpBatchConfiguration.root)

  private val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions

  val monitoringScriptPath: Option[Path] = workflowOptions.get(WorkflowOptionKeys.MonitoringScript).toOption map {
    path =>
      // Fail here if the path exists but can't be built
      getPath(path).get
  }
  override def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): GcpBatchJobPaths =
    new GcpBatchJobPaths(workflowPaths.asInstanceOf[GcpBatchWorkflowPaths], jobKey)
  override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths =
    this.copy(workflowDescriptor = workflowDescriptor)
  override def config: Config = gcpBatchConfiguration.configurationDescriptor.backendConfig
}
