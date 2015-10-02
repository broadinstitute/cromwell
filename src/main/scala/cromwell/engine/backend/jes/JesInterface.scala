package cromwell.engine.backend.jes

import java.net.URL

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.http.HttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.genomics.Genomics
import cromwell.engine.backend.jes.JesInterface.GoogleGenomics
import cromwell.util.google.{GoogleCredentialFactory, GoogleCloudStorage}


object JesInterface {
  val jsonFactory = JacksonFactory.getDefaultInstance
  val httpTransport = GoogleNetHttpTransport.newTrustedTransport

  def apply(appName: String, endpointUrl: URL): JesInterface = {
    val credential = GoogleCredentialFactory.from(jsonFactory, httpTransport)
    val genomics = GoogleGenomics.from(appName, endpointUrl, credential, jsonFactory, httpTransport)
    val storage = GoogleCloudStorage(appName, credential, jsonFactory, httpTransport)

    JesInterface(credential, genomics, storage)
  }

  // Wrapper object around Google's Genomics class providing a convenience 'from' "method"
  object GoogleGenomics {
    def from(applicationName: String, endpointUrl: URL, credential: Credential, jsonFactory: JsonFactory, httpTransport: HttpTransport): Genomics = {
      new Genomics.Builder(httpTransport, jsonFactory, credential).setApplicationName(applicationName).setRootUrl(endpointUrl.toString).build
    }
  }
}

case class JesInterface(credential: Credential, genomics: Genomics, storage: GoogleCloudStorage)