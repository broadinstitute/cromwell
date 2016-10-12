package cromwell.backend.impl.jes

import java.net.URL

import com.google.api.services.genomics.Genomics
import cromwell.core.WorkflowOptions
import cromwell.filesystems.gcs.auth.GoogleAuthMode


case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL) {

  def withOptions(options: WorkflowOptions) = {
    val credential = authMode.credential(options)

    new Genomics.Builder(
      credential.getTransport,
      credential.getJsonFactory,
      credential)
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build
  }
}
