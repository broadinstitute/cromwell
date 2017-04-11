package cromwell.backend.impl.jes

import java.net.URL

import com.google.api.client.http.{HttpRequest, HttpRequestInitializer}
import com.google.api.services.genomics.Genomics
import com.google.auth.http.HttpCredentialsAdapter
import cromwell.core.WorkflowOptions
import cromwell.filesystems.gcs.auth.GoogleAuthMode


case class GenomicsFactory(applicationName: String, authMode: GoogleAuthMode, endpointUrl: URL) {

  def withOptions(options: WorkflowOptions) = {
    val scopedCredentials = authMode.credential(options)
    
    val httpRequestInitializer = {
      val delegate = new HttpCredentialsAdapter(scopedCredentials)
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
