package cromwell.cloudsupport.gcp.genomics

import java.net.URL
import java.time.OffsetDateTime
import java.util.Date

import com.google.api.client.http.GenericUrl
import com.google.auth.oauth2.{AccessToken, GoogleCredentials}
import cromwell.cloudsupport.gcp.auth.MockAuthMode
import org.scalatest.{FlatSpec, Matchers}

class GenomicsFactorySpec extends FlatSpec with Matchers {

  behavior of "GenomicsFactory"

  it should "build a genomics" in {
    val factory = GenomicsFactory("genomics-factory-spec", MockAuthMode, new URL("http://example.com"))
    val accessToken = new AccessToken("my_token", Date.from(OffsetDateTime.now().plusDays(1).toInstant))
    val credentials = GoogleCredentials.of(accessToken)
    val genomics = factory.fromCredentials(credentials)
    genomics.getApplicationName should be("genomics-factory-spec")
    genomics.getRootUrl should be("http://example.com/")
    genomics.getRequestFactory.buildGetRequest(new GenericUrl(new URL("http://example.org")))
  }

}
