package cromwell.filesystems.gcs

import java.nio.file.Paths

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.typesafe.config.ConfigFactory
import cromwell.filesystems.gcs.GoogleAuthMode.EnhancedCredentials
import org.scalatest.{FlatSpec, Matchers}

import scala.util.Try

class GoogleCredentialFactorySpec extends FlatSpec with Matchers {
  import GoogleCredentialFactorySpec._

  behavior of "GoogleCredentialFactory"

  it should "refresh a token using user credentials" taggedAs GcsIntegrationTest in {
    val credential = UserMode(
      name = "user",
      user = secretConf("user"),
      secretsFile = secretConf("secrets-file"),
      datastoreDir = secretConf("data-store-dir")).credential(emptyOptions)

    val firstCredentialTry: Try[Credential] = credential.freshCredential
    assert(firstCredentialTry.isSuccess)
    val firstCredential = firstCredentialTry.get
    firstCredential.getAccessToken shouldNot be(empty)

    firstCredential.setExpiresInSeconds(59L)

    val secondCredentialTry: Try[Credential] = firstCredential.freshCredential
    assert(secondCredentialTry.isSuccess)

    val secondCredential = secondCredentialTry.get
    secondCredential.getAccessToken shouldNot be(empty)
    secondCredential.getExpiresInSeconds shouldNot be(null)
    secondCredential.getExpiresInSeconds.longValue should be > 60L
  }

  it should "refresh a token using a service account" taggedAs GcsIntegrationTest in {
    val credential = ServiceAccountMode(
      name = "service",
      accountId = secretConf("service-account-id"),
      pemPath = secretConf("pem-file")).credential(emptyOptions)

    val firstCredentialTry: Try[Credential] = credential.freshCredential
    assert(firstCredentialTry.isSuccess)
    val firstCredential = firstCredentialTry.get
    firstCredential.getAccessToken shouldNot be(empty)

    firstCredential.setExpiresInSeconds(59L)

    val secondCredentialTry: Try[Credential] = firstCredential.freshCredential
    assert(secondCredentialTry.isSuccess)

    val secondCredential = secondCredentialTry.get
    secondCredential.getAccessToken shouldNot be(empty)
    secondCredential.getExpiresInSeconds shouldNot be(null)
    secondCredential.getExpiresInSeconds.longValue should be > 60L
  }

  it should "refresh a token using a refresh token" taggedAs GcsIntegrationTest in {
    val opts = GoogleOptionsMap(Map("refresh_token" -> secretConf("refresh_token")))

    val credential = RefreshTokenMode(name = "refresh",
      clientId = secretConf("client-id"),
      clientSecret = secretConf("client-secret")).credential(opts)

    val firstUserCredentialsTry = credential.freshCredential

    assert(firstUserCredentialsTry.isSuccess)
    val firstUserCredentials = firstUserCredentialsTry.get

    val firstRefreshedUserCredentialsTry: Try[Credential] = firstUserCredentials.freshCredential
    assert(firstRefreshedUserCredentialsTry.isSuccess)
    val firstRefreshedUserCredentials = firstRefreshedUserCredentialsTry.get
    firstRefreshedUserCredentials.getAccessToken shouldNot be(empty)

    firstRefreshedUserCredentials.setExpiresInSeconds(59L)

    val secondRefreshedUserCredentialsTry: Try[Credential] = firstRefreshedUserCredentials.freshCredential
    assert(secondRefreshedUserCredentialsTry.isSuccess)

    val secondRefreshedUserCredentials = secondRefreshedUserCredentialsTry.get
    secondRefreshedUserCredentials.getAccessToken shouldNot be(empty)
    secondRefreshedUserCredentials.getExpiresInSeconds shouldNot be(null)
    secondRefreshedUserCredentials.getExpiresInSeconds.longValue should be > 60L
  }

  it should "not refresh an empty token" in {

    val wrongCredentials = new GoogleCredential.Builder()
      .setTransport(GoogleNetHttpTransport.newTrustedTransport)
      .setJsonFactory(JacksonFactory.getDefaultInstance)
      .setClientSecrets("fakeId", "fakeSecret")
      .build()

    val exception = wrongCredentials.freshCredential.failed.get

    exception.getMessage should be("Unable to refresh token")
  }

  it should "refresh a token using application default credentials" taggedAs GcsIntegrationTest in {
    val credential = applicationDefaultCredential

    val firstCredentialTry: Try[Credential] = credential.freshCredential
    assert(firstCredentialTry.isSuccess)
    val firstCredential = firstCredentialTry.get
    firstCredential.getAccessToken shouldNot be(empty)

    firstCredential.setExpiresInSeconds(59L)

    val secondCredentialTry: Try[Credential] = firstCredential.freshCredential
    assert(secondCredentialTry.isSuccess)

    val secondCredential = secondCredentialTry.get
    secondCredential.getAccessToken shouldNot be(empty)
    secondCredential.getExpiresInSeconds shouldNot be(null)
    secondCredential.getExpiresInSeconds.longValue should be > 60L
  }
}

object GoogleCredentialFactorySpec {
  /*

  To run this integration spec, your cromwell-credentials.conf file should have the following keys for the listed tests:

  // For testing UserMode
  user = "<user email>"
  secrets-file = "<path to json secrets for above email>"
  data-store-dir = "<path to where secrets will be stored>"

  // For testing ServiceAccountMode
  service-account-id = "<service account id>"
  pem-file = "<path to the pem file>"

  // For testing RefreshTokenMode
  client-id = "<service account id, preferably one configured to use the oauth playground>"
  client-secret = "<secret for the above client-id>"
  refresh_token = "<refresh token created with the above id/secret and the cloud storage api scopes>"

   */

  private lazy val credentialsConfig = ConfigFactory.parseFile(Paths.get("cromwell-credentials.conf").toFile)

  private def secretConf(path: String) = credentialsConfig.getString(path)

  private val emptyOptions = GoogleOptionsMap(Map.empty)

  def applicationDefaultCredential = ApplicationDefaultMode(name = "default").credential(emptyOptions)
}

case class GoogleOptionsMap(map: Map[String, String]) extends GoogleAuthMode.GoogleAuthOptions {
  override def get(key: String): Try[String] = Try { map(key) }
}
