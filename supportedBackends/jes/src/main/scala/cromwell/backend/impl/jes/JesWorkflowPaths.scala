package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.backend.impl.jes.authentication.JesCredentials
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.WorkflowOptions.FinalCallLogsDir
import cromwell.filesystems.gcs.{GcsFileSystem, GcsFileSystemProvider, GoogleAuthMode}

import scala.concurrent.ExecutionContext

object JesWorkflowPaths {
  private val GcsRootOptionKey = "jes_gcs_root"
  private val AuthFilePathOptionKey = "auth_bucket"

  def apply(workflowDescriptor: BackendWorkflowDescriptor,
            jesConfiguration: JesConfiguration,
            credentials: JesCredentials)(implicit ec: ExecutionContext) = {
    new JesWorkflowPaths(workflowDescriptor, jesConfiguration, credentials)
  }
}

class JesWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor,
                       jesConfiguration: JesConfiguration,
                       credentials: JesCredentials)(implicit ec: ExecutionContext) {

  private val gcsStorage = GoogleAuthMode.buildStorage(credentials.gcsCredential, jesConfiguration.googleConfig.applicationName)
  val gcsFileSystemProvider: GcsFileSystemProvider = GcsFileSystemProvider(gcsStorage)(ec)
  val gcsFileSystem = GcsFileSystem(gcsFileSystemProvider)

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
    val genomicsStorage = GoogleAuthMode.buildStorage(credentials.genomicsCredential, jesConfiguration.googleConfig.applicationName)
    val fileSystemWithGenomicsAuth = GcsFileSystem(GcsFileSystemProvider(genomicsStorage)(ec))
    val bucket = workflowDescriptor.workflowOptions.get(JesWorkflowPaths.AuthFilePathOptionKey) getOrElse workflowRootPath.toString

    fileSystemWithGenomicsAuth.getPath(bucket).resolve(s"${workflowDescriptor.id}_auth.json")
  }

  def toJesCallPaths(jobKey: BackendJobDescriptorKey) = JesCallPaths(jobKey, workflowDescriptor, jesConfiguration, credentials)(ec)
}
