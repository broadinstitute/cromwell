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
import PipelinesApiWorkflowPaths._

import scala.concurrent.ExecutionContext
import scala.language.postfixOps

object PipelinesApiWorkflowPaths {
  private val GcsRootOptionKey = "jes_gcs_root"
  private val AuthFilePathOptionKey = "auth_bucket"
  private val GcsPrefix = "gs://"

  private[common] def callCachePathPrefixFromExecutionRoot(executionRoot: String): String = {
    // If the root looks like gs://bucket/stuff-under-bucket this should return gs://bucket
    GcsPrefix + executionRoot.substring(GcsPrefix.length).takeWhile(_ != '/')
  }
}

case class PipelinesApiWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor,
                                     gcsCredentials: Credentials,
                                     genomicsCredentials: Credentials,
                                     papiConfiguration: PipelinesApiConfiguration,
                                     override val pathBuilders: PathBuilders,
                                     // This allows for the adjustment of the standard stream file names in PAPI v1 to match the
                                     // combined controller + job standard output and error files. PAPI v1 controls the periodic
                                     // delocalization of these files so the metadata Cromwell publishes for these files needs
                                     // to match the PAPI v1 names.
                                     standardStreamNameToFileNameMetadataMapper: (PipelinesApiJobPaths, String) => String)(implicit ec: ExecutionContext) extends WorkflowPaths {

  override lazy val executionRootString: String =
    workflowDescriptor.workflowOptions.getOrElse(PipelinesApiWorkflowPaths.GcsRootOptionKey, papiConfiguration.root)

  override lazy val callCacheRootPrefix: Option[String] = Option(callCachePathPrefixFromExecutionRoot(executionRootString))

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
      Option(papiConfiguration.papiAttributes.project)
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
