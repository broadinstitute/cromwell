package cromwell.backend.impl.jes

import akka.actor.ActorSystem
import com.google.auth.Credentials
import com.typesafe.config.Config
import cromwell.backend.impl.jes.JesAsyncBackendJobExecutionActor.WorkflowOptionKeys
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.WorkflowOptions
import cromwell.core.path.{Path, PathBuilder}
import cromwell.filesystems.gcs.{GcsPathBuilder, GcsPathBuilderFactory}

import scala.language.postfixOps

object JesWorkflowPaths {
  private val GcsRootOptionKey = "jes_gcs_root"
  private val AuthFilePathOptionKey = "auth_bucket"
}

case class JesWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor,
                            gcsCredentials: Credentials,
                            genomicsCredentials: Credentials,
                            jesConfiguration: JesConfiguration)(implicit actorSystem: ActorSystem) extends WorkflowPaths {

  override lazy val executionRootString: String =
    workflowDescriptor.workflowOptions.getOrElse(JesWorkflowPaths.GcsRootOptionKey, jesConfiguration.root)

  private val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions

  private val gcsPathBuilder: GcsPathBuilder = jesConfiguration.gcsPathBuilderFactory.fromCredentials(workflowOptions, gcsCredentials)

  val gcsAuthFilePath: Path = {
    // The default auth file bucket is always at the root of the root workflow
    val defaultBucket = executionRoot.resolve(workflowDescriptor.rootWorkflow.unqualifiedName).resolve(workflowDescriptor.rootWorkflowId.toString)
    val bucket = workflowDescriptor.workflowOptions.get(JesWorkflowPaths.AuthFilePathOptionKey) getOrElse defaultBucket.pathAsString

    /*
     * This is an "exception". The filesystem used here is built from genomicsAuth
     * unlike everywhere else where the filesystem used is built from gcsFileSystemAuth
     */
    val pathBuilderWithGenomicsAuth = GcsPathBuilder.fromCredentials(
      genomicsCredentials,
      jesConfiguration.googleConfig.applicationName,
      None,
      GcsPathBuilderFactory.DefaultCloudStorageConfiguration,
      workflowOptions
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

  override def toJobPaths(workflowPaths: WorkflowPaths, jobKey: BackendJobDescriptorKey): JesJobPaths = {
    new JesJobPaths(workflowPaths.asInstanceOf[JesWorkflowPaths], jobKey)
  }

  override protected def withDescriptor(workflowDescriptor: BackendWorkflowDescriptor): WorkflowPaths = this.copy(workflowDescriptor = workflowDescriptor)

  override def config: Config = jesConfiguration.configurationDescriptor.backendConfig
  override def pathBuilders: List[PathBuilder] = List(gcsPathBuilder)
}
