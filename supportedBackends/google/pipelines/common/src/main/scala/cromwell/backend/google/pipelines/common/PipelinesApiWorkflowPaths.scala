package cromwell.backend.google.pipelines.common

import com.google.api.gax.retrying.RetrySettings
import com.google.auth.Credentials
import com.typesafe.config.Config
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.cloudsupport.gcp.gcs.GcsStorage
import cromwell.core.WorkflowOptions
import cromwell.core.path.Path
import cromwell.core.path.PathFactory.PathBuilders
import cromwell.filesystems.gcs.GcsPathBuilder
import cromwell.filesystems.gcs.cache.GcsBucketInformationPolicies

import scala.concurrent.ExecutionContext
import scala.language.postfixOps

object PipelinesApiWorkflowPaths {
  private val GcsRootOptionKey = "jes_gcs_root"
  private val AuthFilePathOptionKey = "auth_bucket"
}

case class PipelinesApiWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor,
                                     gcsCredentials: Credentials,
                                     genomicsCredentials: Credentials,
                                     papiConfiguration: PipelinesApiConfiguration,
                                     override val pathBuilders: PathBuilders)(implicit ec: ExecutionContext) extends WorkflowPaths {
  override lazy val executionRootString: String =
    workflowDescriptor.workflowOptions.getOrElse(PipelinesApiWorkflowPaths.GcsRootOptionKey, papiConfiguration.root)

  private val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions

  val gcsAuthFilePath: Path = {
    // The default auth file bucket is always at the root of the root workflow
    val defaultBucket = executionRoot.resolve(workflowDescriptor.rootWorkflow.name).resolve(workflowDescriptor.rootWorkflowId.toString)
    val bucket = workflowDescriptor.workflowOptions.get(PipelinesApiWorkflowPaths.AuthFilePathOptionKey) getOrElse defaultBucket.pathAsString

    /*
     * This is an "exception". The filesystem used here is built from genomicsAuth
     * unlike everywhere else where the filesystem used is built from gcsFileSystemAuth
     */
    val pathBuilderWithGenomicsAuth = GcsPathBuilder.fromCredentials(
      genomicsCredentials,
      papiConfiguration.googleConfig.applicationName,
      RetrySettings.newBuilder().build(),
      GcsStorage.DefaultCloudStorageConfiguration,
      workflowOptions,
      Option(papiConfiguration.jesAttributes.project),
      // This is only used for Pipelines V1, no need to go out of our way to support requester pays there, so use the static default
      GcsBucketInformationPolicies.OnDemandPolicy
    )

    val authBucket = pathBuilderWithGenomicsAuth.build(bucket) recover {
      case ex => throw new Exception(s"Invalid gcs auth_bucket path $bucket", ex)
    } get

    authBucket.resolve(s"${workflowDescriptor.rootWorkflowId}_auth.json")
  }

  val monitoringScriptPath: Option[Path] = workflowOptions.get(WorkflowOptionKeys.MonitoringScript).toOption map { path =>
    // Fail here if the path exists but can't be built
    getPath(path).get
  }

  override def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): PipelinesApiJobPaths = {
    new PipelinesApiJobPaths(workflowPaths.asInstanceOf[PipelinesApiWorkflowPaths], jobKey)
  }

  override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = this.copy(workflowDescriptor = workflowDescriptor)

  override def config: Config = papiConfiguration.configurationDescriptor.backendConfig
}
