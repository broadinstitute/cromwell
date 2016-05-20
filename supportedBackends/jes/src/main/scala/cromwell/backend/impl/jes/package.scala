package cromwell.backend.impl

import cromwell.backend.impl.jes.JesImplicits.GoogleAuthWorkflowOptions
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.filesystems.gcs.{GcsFileSystem, GcsFileSystemProvider, GoogleAuthMode, GoogleConfiguration}

package object jes {
  def buildGcsFileSystem(configurationDescriptor: BackendConfigurationDescriptor, workflowDescriptor: BackendWorkflowDescriptor): GcsFileSystem = {
      // PBE what follows lacks any semblance of error checking

      // PBE check the config sanity in the actory factory, preferably in a fail-fast manner.
      val genomicsAuth: GoogleAuthMode =
        GoogleConfiguration(configurationDescriptor.globalConfig).auth(configurationDescriptor.backendConfig.getString("genomics.auth")) getOrElse { throw new RuntimeException("borked config") }

      // PBE It might be nice to only build this filesystems once per workflow (maybe in the initialization actor) and
      // somehow make that available to the workflow's job execution actors.
      val authOptions = workflowDescriptor.workflowOptions.toGoogleAuthOptions
      GcsFileSystemProvider(genomicsAuth.buildStorage(authOptions, configurationDescriptor.globalConfig)).getFileSystem
  }
}
