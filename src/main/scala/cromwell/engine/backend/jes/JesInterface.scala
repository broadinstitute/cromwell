package cromwell.engine.backend.jes

import java.net.URL

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.http.HttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.client.json.jackson2.JacksonFactory
import com.google.api.services.genomics.Genomics
import cromwell.util.google.{GoogleCloudStorage, GoogleCredentialFactory}


object JesInterface {
  def apply(appName: String): JesInterface = {
    val jsonFactory = JacksonFactory.getDefaultInstance
    val httpTransport = GoogleNetHttpTransport.newTrustedTransport
    val credential = GoogleCredentialFactory.from(jsonFactory, httpTransport)
    val genomics = GoogleGenomics.from(appName, credential, jsonFactory, httpTransport)
    val storage = GoogleCloudStorage(appName, credential, jsonFactory, httpTransport)

    JesInterface(credential, genomics, storage)
  }

  // Wrapper object around Google's Genomics class providing a convenience 'from' "method"
  object GoogleGenomics {
    def GenomicsUrl = new URL("https://staging-genomics.sandbox.googleapis.com")

    def from(applicationName: String, credential: Credential, jsonFactory: JsonFactory, httpTransport: HttpTransport): Genomics = {
      new Genomics.Builder(httpTransport, jsonFactory, credential).setApplicationName(applicationName).setRootUrl(GenomicsUrl.toString).build
    }
  }
}

case class JesInterface(credential: Credential, genomics: Genomics, storage: GoogleCloudStorage)
