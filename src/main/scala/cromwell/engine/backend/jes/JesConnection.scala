package cromwell.engine.backend.jes

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.http.HttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.services.genomics.Genomics
import cromwell.util.google.{GoogleCloudStorage, GoogleCredential, GoogleGenomics}


object JesConnection {
  def apply(appName: String, jsonFactory: JsonFactory, httpTransport: HttpTransport): JesConnection = {
    val credential = GoogleCredential.from(jsonFactory, httpTransport)
    val genomics = GoogleGenomics.from(appName, credential, jsonFactory, httpTransport)
    val storage = GoogleCloudStorage(appName, credential, jsonFactory, httpTransport)

    JesConnection(credential, genomics, storage)
  }
}

case class JesConnection(credential: Credential, genomics: Genomics, storage: GoogleCloudStorage)
