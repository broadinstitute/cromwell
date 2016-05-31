package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.backend.impl.jes.io._

object JesWorkflowPaths {
  private val GcsRootOptionKey = "jes_gcs_root"
  private val AuthFilePathOptionKey = "auth_bucket"

  def apply(workflowDescriptor: BackendWorkflowDescriptor, jesConfiguration: JesConfiguration) = {
    new JesWorkflowPaths(workflowDescriptor, jesConfiguration)
  }
}

class JesWorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor, jesConfiguration: JesConfiguration) {
  val gcsFileSystem = buildFilesystem(workflowDescriptor, jesConfiguration.jesAttributes.gcsFilesystemAuth, jesConfiguration.googleConfig)

  val rootPath: Path =
    gcsFileSystem.getPath(workflowDescriptor.workflowOptions.getOrElse(JesWorkflowPaths.GcsRootOptionKey, jesConfiguration.root))

  val workflowRootPath: Path = rootPath.resolve(workflowDescriptor.workflowNamespace.workflow.unqualifiedName)
    .resolve(workflowDescriptor.id.toString)

  val gcsAuthFilePath: Path = {
    /*
     * This is an "exception". The filesystem used here is built from genomicsAuth
     * unlike everywhere else where the filesystem used is built from gcsFileSystemAuth
     */
    val fs = buildFilesystem(workflowDescriptor, jesConfiguration.jesAttributes.genomicsAuth, jesConfiguration.googleConfig)
    val bucket = workflowDescriptor.workflowOptions.get(JesWorkflowPaths.AuthFilePathOptionKey) getOrElse workflowRootPath.toString
    fs.getPath(bucket).resolve(s"${workflowDescriptor.id}_auth.json")
  }
}
