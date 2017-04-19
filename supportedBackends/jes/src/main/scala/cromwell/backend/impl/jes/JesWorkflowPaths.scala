package cromwell.backend.impl.jes

import akka.actor.ActorSystem
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

  def apply(workflowDescriptor: BackendWorkflowDescriptor,
            jesConfiguration: JesConfiguration)(implicit actorSystem: ActorSystem): JesWorkflowPaths = {
    new JesWorkflowPaths(workflowDescriptor, jesConfiguration)
  }
}

class JesWorkflowPaths(val workflowDescriptor: BackendWorkflowDescriptor,
                       jesConfiguration: JesConfiguration)(implicit actorSystem: ActorSystem) extends WorkflowPaths {

  override lazy val executionRootString: String =
    workflowDescriptor.workflowOptions.getOrElse(JesWorkflowPaths.GcsRootOptionKey, jesConfiguration.root)
  private val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions
  val gcsPathBuilder: GcsPathBuilder = jesConfiguration.gcsPathBuilderFactory.withOptions(workflowOptions)

  val gcsAuthFilePath: Path = {
    /*
     * This is an "exception". The filesystem used here is built from genomicsAuth
     * unlike everywhere else where the filesystem used is built from gcsFileSystemAuth
     */
    val genomicsCredentials = jesConfiguration.jesAuths.genomics

    // The default auth file bucket is always at the root of the root workflow
    val defaultBucket = executionRoot.resolve(workflowDescriptor.rootWorkflow.unqualifiedName).resolve(workflowDescriptor.rootWorkflowId.toString)

    val bucket = workflowDescriptor.workflowOptions.get(JesWorkflowPaths.AuthFilePathOptionKey) getOrElse defaultBucket.pathAsString
    val authBucket = GcsPathBuilderFactory(genomicsCredentials, jesConfiguration.googleConfig.applicationName).withOptions(workflowOptions).build(bucket) recover {
      case ex => throw new Exception(s"Invalid gcs auth_bucket path $bucket", ex)
    } get

    authBucket.resolve(s"${workflowDescriptor.rootWorkflowId}_auth.json")
  }

  val monitoringPath: Option[Path] = workflowOptions.get(WorkflowOptionKeys.MonitoringScript).toOption map { path =>
    // Fail here if the path exists but can't be built
    getPath(path).get
  }

  override def toJobPaths(jobKey: BackendJobDescriptorKey, jobWorkflowDescriptor: BackendWorkflowDescriptor): JesJobPaths = JesJobPaths(this, jobKey)
  override def config: Config = jesConfiguration.configurationDescriptor.backendConfig
  override def pathBuilders: List[PathBuilder] = List(gcsPathBuilder)
}
