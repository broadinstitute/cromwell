package cromwell.engine.backend.jes

import java.net.URL

import com.google.api.services.genomics.Genomics
import cromwell.filesystems.gcs.{GoogleCredentialFactory, GoogleConfiguration}
import cromwell.filesystems.gcs.GoogleCredentialFactory
import cromwell.util.google.GenomicsScopes

object GenomicsFactory {
  def apply(genomicsConf: GoogleConfiguration, endpointUrl: URL): Genomics = {
    val credential = GoogleCredentialFactory(genomicsConf.authMode, GenomicsScopes.genomicsScopes)
    new Genomics.Builder(credential.getTransport, credential.getJsonFactory, credential)
      .setApplicationName(genomicsConf.appName)
      .setRootUrl(endpointUrl.toString)
      .build
  }
}
