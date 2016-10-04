package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.ActorSystem
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.WorkflowOptions
import cromwell.core.WorkflowOptions.FinalCallLogsDir
import cromwell.filesystems.gcs.{RetryableGcsPathBuilder, GcsPathBuilder, GcsPathBuilderFactory}

import scala.language.postfixOps
import scala.util.Try

object JesWorkflowPaths {
  private val GcsRootOptionKey = "jes_gcs_root"
  private val AuthFilePathOptionKey = "auth_bucket"

  def apply(workflowDescriptor: BackendWorkflowDescriptor,
            jesConfiguration: JesConfiguration)(implicit actorSystem: ActorSystem) = {
    new JesWorkflowPaths(workflowDescriptor, jesConfiguration)
  }
}

class JesWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor,
                       jesConfiguration: JesConfiguration)(implicit actorSystem: ActorSystem) {

  private val rootString = workflowDescriptor.workflowOptions.getOrElse(JesWorkflowPaths.GcsRootOptionKey, jesConfiguration.root)
  private val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions

  val gcsPathBuilder: RetryableGcsPathBuilder = jesConfiguration.gcsPathBuilderFactory.withOptions(workflowOptions)

  def getPath(gcsUrl: String): Try[CloudStoragePath] = gcsPathBuilder.build(gcsUrl)
  def getHash(gcsUrl: CloudStoragePath) = gcsPathBuilder.getHash(gcsUrl)

  val rootPath: Path = getPath(rootString) recover {
    case ex => throw new Exception(s"Failed to : $rootString", ex)
  } get

  val workflowRootPath: Path = rootPath.resolve(workflowDescriptor.workflowNamespace.workflow.unqualifiedName).resolve(s"${workflowDescriptor.id.toString}/")

  val finalCallLogsPath = workflowDescriptor.getWorkflowOption(FinalCallLogsDir) map getPath map { _.get }

  val gcsAuthFilePath: Path = {
    /*
     * This is an "exception". The filesystem used here is built from genomicsAuth
     * unlike everywhere else where the filesystem used is built from gcsFileSystemAuth
     */
    val genomicsCredentials = jesConfiguration.jesAuths.genomics
    val bucket = workflowDescriptor.workflowOptions.get(JesWorkflowPaths.AuthFilePathOptionKey) getOrElse workflowRootPath.toUri.toString
    val authBucket = GcsPathBuilderFactory(genomicsCredentials).withOptions(workflowOptions).build(bucket) recover {
      case ex => throw new Exception(s"Invalid gcs auth_bucket path $bucket", ex)
    } get

    authBucket.resolve(s"${workflowDescriptor.id}_auth.json")
  }

  def toJesCallPaths(jobKey: BackendJobDescriptorKey) = JesCallPaths(jobKey, workflowDescriptor, jesConfiguration)
}
