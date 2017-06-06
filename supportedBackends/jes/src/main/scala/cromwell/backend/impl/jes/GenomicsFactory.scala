package cromwell.backend.impl.jes

import java.net.URL

import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.services.genomics.Genomics
import com.google.auth.Credentials
import com.google.auth.http.HttpCredentialsAdapter
import cromwell.filesystems.gcs.auth.GoogleAuthMode


case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL) {

  def fromCredentials(credentials: Credentials) = {
    val httpRequestInitializer = {
      val delegate = new HttpCredentialsAdapter(credentials)
      new HttpRequestInitializer() {
        def initialize(httpRequest: HttpRequest) = {
          delegate.initialize(httpRequest)
        }
      }
    }

    new Genomics.Builder(
      GoogleAuthMode.httpTransport,
      GoogleAuthMode.jsonFactory,
      httpRequestInitializer)
      .setApplicationName(applicationName)
      .setRootUrl(endpointUrl.toString)
      .build
  }
}
