package cromwell.util.google

import java.net.URL

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.http.HttpTransport
import com.google.api.client.json.JsonFactory
import com.google.api.services.genomics.Genomics

object GoogleGenomics {
  def GenomicsUrl = new URL("https://staging-genomics.sandbox.googleapis.com")

  def from(applicationName: String, credential: Credential, jsonFactory: JsonFactory, httpTransport: HttpTransport): Genomics = {
    new Genomics.Builder(httpTransport, jsonFactory, credential).setApplicationName(applicationName).setRootUrl(GenomicsUrl.toString).build
  }
}
