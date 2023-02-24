package cromwell.backend.google.pipelines.batch

//import com.google.auth.Credentials
import com.typesafe.config.Config
import cromwell.backend.google.pipelines.batch.GcpBatchWorkflowPaths.callCachePathPrefixFromExecutionRoot
import cromwell.backend.google.pipelines.common.WorkflowOptionKeys
import cromwell.backend.io.WorkflowPaths
import cromwell.core.WorkflowOptions
import cromwell.core.path.Path

//import cromwell.core.WorkflowOptions

//import scala.concurrent.ExecutionContext
//import cromwell.backend.io.{JobPaths, WorkflowPaths}
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
//import cromwell.backend.google.pipelines.common.PipelinesApiWorkflowPaths
//import cromwell.backend.io.WorkflowPaths
//import cromwell.core.WorkflowOptions
import cromwell.core.path.PathFactory.PathBuilders


object GcpBatchWorkflowPaths {
  private val GcsRootOptionKey = "gcp_batch_gcs_root"
  //private val AuthFilePathOptionKey = "auth_bucket"
  private val GcsPrefix = "gs://"

  private def callCachePathPrefixFromExecutionRoot(executionRoot: String): String = {
    // If the root looks like gs://bucket/stuff-under-bucket this should return gs://bucket
    GcsPrefix + executionRoot.substring(GcsPrefix.length).takeWhile(_ != '/')
  }
}
case class GcpBatchWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor,
                                 gcpBatchConfiguration: GcpBatchConfiguration,
                                 override val pathBuilders: PathBuilders) extends WorkflowPaths {

  override lazy val executionRootString: String = workflowDescriptor.workflowOptions.getOrElse(GcpBatchWorkflowPaths.GcsRootOptionKey, gcpBatchConfiguration.root)
  override lazy val callCacheRootPrefix: Option[String] = Option(callCachePathPrefixFromExecutionRoot(executionRootString))

  private val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions

  val monitoringScriptPath: Option[Path] = workflowOptions.get(WorkflowOptionKeys.MonitoringScript)
                                                          .toOption map { path =>
    // Fail here if the path exists but can't be built
    getPath(path).get
  }
  override def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): GcpBatchJobPaths = {
    new GcpBatchJobPaths(workflowPaths.asInstanceOf[GcpBatchWorkflowPaths], jobKey)
  }
  override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = this.copy(workflowDescriptor = workflowDescriptor)
  def config: Config = gcpBatchConfiguration.configurationDescriptor.backendConfig
}



