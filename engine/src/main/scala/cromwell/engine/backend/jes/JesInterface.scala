package cromwell.engine.backend.jes

import java.net.URL

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.http.HttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.services.genomics.Genomics
import cromwell.engine.backend.io.filesystem.gcs.GcsFileSystem
import cromwell.engine.io.gcs.GoogleConfiguration
import cromwell.util.google.GoogleCredentialFactory

case class JesInterface(gcsFileSystem: GcsFileSystem, genomics: Genomics)


object GenomicsFactory {

  def apply(googleConfiguration: GoogleConfiguration, endpointUrl: URL): Genomics = {
    val credential = GoogleCredentialFactory.fromCromwellAuthScheme(googleConfiguration)
    GoogleGenomics.from(googleConfiguration.appName, endpointUrl, credential, credential.getJsonFactory, credential.getTransport)
  }

  // Wrapper object around Google's Genomics class providing a convenience 'from' "method"
  object GoogleGenomics {
    def from(applicationName: String, endpointUrl: URL, credential: Credential, jsonFactory: JsonFactory, httpTransport: HttpTransport): Genomics = {
      new Genomics.Builder(httpTransport, jsonFactory, credential).setApplicationName(applicationName).setRootUrl(endpointUrl.toString).build
    }
  }
}
