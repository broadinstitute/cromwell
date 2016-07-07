package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.WorkflowOptions
import cromwell.core.WorkflowOptions.FinalCallLogsDir
import cromwell.filesystems.gcs.GoogleAuthMode.GoogleAuthOptions
import cromwell.filesystems.gcs.{GcsFileSystem, GcsFileSystemProvider, GoogleAuthMode}

import scala.util.Try

object JesWorkflowPaths {
  private val GcsRootOptionKey = "jes_gcs_root"
  private val AuthFilePathOptionKey = "auth_bucket"

  def apply(workflowDescriptor: BackendWorkflowDescriptor, jesConfiguration: JesConfiguration, gcsFileSystem: GcsFileSystem) = {
    new JesWorkflowPaths(workflowDescriptor, jesConfiguration, gcsFileSystem)
  }

  implicit class GoogleAuthWorkflowOptions(val workflowOptions: WorkflowOptions) extends AnyVal {
    def toGoogleAuthOptions: GoogleAuthMode.GoogleAuthOptions = new GoogleAuthOptions {
      override def get(key: String): Try[String] = workflowOptions.get(key)
    }
  }
}

class JesWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor, jesConfiguration: JesConfiguration, val gcsFileSystem: GcsFileSystem) {
import JesWorkflowPaths._
  val rootPath: Path =
    gcsFileSystem.getPath(workflowDescriptor.workflowOptions.getOrElse(JesWorkflowPaths.GcsRootOptionKey, jesConfiguration.root))

  val workflowRootPath: Path = rootPath.resolve(workflowDescriptor.workflowNamespace.workflow.unqualifiedName)
    .resolve(workflowDescriptor.id.toString)

  val finalCallLogsPath = workflowDescriptor.getWorkflowOption(FinalCallLogsDir) map { gcsFileSystem.getPath(_) }

  val gcsAuthFilePath: Path = {
    /*
     * This is an "exception". The filesystem used here is built from genomicsAuth
     * unlike everywhere else where the filesystem used is built from gcsFileSystemAuth
     */
    val bucket = workflowDescriptor.workflowOptions.get(JesWorkflowPaths.AuthFilePathOptionKey) getOrElse workflowRootPath.toString

    val storage = jesConfiguration.jesAttributes.genomicsAuth.buildStorage(workflowDescriptor.workflowOptions.toGoogleAuthOptions, jesConfiguration.googleConfig)
    val fileSystemWithGenomicsAuth = GcsFileSystem(GcsFileSystemProvider(storage))

    fileSystemWithGenomicsAuth.getPath(bucket).resolve(s"${workflowDescriptor.id}_auth.json")
  }

  def toJesCallPaths(jobKey: BackendJobDescriptorKey) = JesCallPaths(jobKey, workflowDescriptor, jesConfiguration, gcsFileSystem)
}
