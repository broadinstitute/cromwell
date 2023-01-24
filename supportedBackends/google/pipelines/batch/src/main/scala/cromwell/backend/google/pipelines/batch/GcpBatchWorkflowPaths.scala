package cromwell.backend.google.pipelines.batch

//import com.google.auth.Credentials
import com.typesafe.config.Config
import cromwell.backend.io.WorkflowPaths
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
  //private val GcsPrefix = "gs://"
}
//case class GcpBatchWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor, gcsCredentials: Credentials, gcpBatchConfiguration: GcpBatchConfiguration, override val pathBuilders: PathBuilders) extends WorkflowPaths {
case class GcpBatchWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor,
                                 gcpBatchConfiguration: GcpBatchConfiguration,
                                 override val pathBuilders: PathBuilders) extends WorkflowPaths {

  override lazy val executionRootString: String = workflowDescriptor.workflowOptions.getOrElse(GcpBatchWorkflowPaths.GcsRootOptionKey, gcpBatchConfiguration.root)
  //override lazy val callCacheRootPrefix: Option[String] = Option(callCachePathPrefixFromExecutionRoot(executionRootString))

  //private val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions

  /*
  val gcsAuthFilePath: Path = {
    // The default auth file bucket is always at the root of the root workflow
    val defaultBucket = executionRoot
      .resolve(workflowDescriptor
        .rootWorkflow
        .name)
      .resolve(workflowDescriptor
        .rootWorkflowId
        .toString)
    val bucket = workflowDescriptor
      .workflowOptions
      .get(GcpBatchWorkflowPaths
        .AuthFilePathOptionKey) getOrElse defaultBucket
      .pathAsString
*/
  override def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): GcpBatchJobPaths = {
    new GcpBatchJobPaths(workflowPaths.asInstanceOf[GcpBatchWorkflowPaths], jobKey)
  }
  //override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = this.copy(workflowDescriptor = workflowDescriptor)

  def config: Config = gcpBatchConfiguration.configurationDescriptor.backendConfig



  override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = ???
}



