package cromwell.engine.backend.jes

import java.net.URL

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.http.HttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.services.genomics.Genomics
import cromwell.util.google.{GoogleCloudStorage, GoogleCredentialFactory}

case class JesInterface(storage: GoogleCloudStorage, genomics: Genomics)

object GcsFactory {
  def apply(appName: String, credential: Credential): GoogleCloudStorage = {
    GoogleCloudStorage(appName, credential, GoogleCredentialFactory.jsonFactory, GoogleCredentialFactory.httpTransport)
  }
}

object GenomicsFactory {
  def apply(appName: String, endpointUrl: URL, credential: Credential): Genomics = {
    GoogleGenomics.from(appName, endpointUrl, credential, GoogleCredentialFactory.jsonFactory, GoogleCredentialFactory.httpTransport)
  }

  // Wrapper object around Google's Genomics class providing a convenience 'from' "method"
  object GoogleGenomics {
    def from(applicationName: String, endpointUrl: URL, credential: Credential, jsonFactory: JsonFactory, httpTransport: HttpTransport): Genomics = {
      new Genomics.Builder(httpTransport, jsonFactory, credential).setApplicationName(applicationName).setRootUrl(endpointUrl.toString).build
    }
  }
}
