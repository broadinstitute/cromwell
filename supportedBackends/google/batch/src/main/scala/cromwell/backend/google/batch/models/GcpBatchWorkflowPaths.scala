package cromwell.backend.google.batch.models

import com.google.auth.Credentials
import com.typesafe.config.Config
import cromwell.backend.google.batch.models.GcpBatchWorkflowPaths.callCachePathPrefixFromExecutionRoot
import cromwell.backend.google.batch.runnable.WorkflowOptionKeys
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.WorkflowOptions
import cromwell.core.path.Path
import cromwell.core.path.PathFactory.PathBuilders

import scala.concurrent.ExecutionContext

object GcpBatchWorkflowPaths {
  val GcsRootOptionKey = "gcp_batch_gcs_root"
  private val GcsPrefix = "gs://"

  private[models] def callCachePathPrefixFromExecutionRoot(executionRoot: String): String =
    // If the root looks like gs://bucket/stuff-under-bucket this should return gs://bucket
    GcsPrefix + executionRoot.substring(GcsPrefix.length).takeWhile(_ != '/')
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
  override lazy val callCacheRootPrefix: Option[String] = Option(
    callCachePathPrefixFromExecutionRoot(executionRootString)
  )

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
