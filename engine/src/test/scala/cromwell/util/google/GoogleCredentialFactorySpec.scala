package cromwell.util.google

import java.nio.file.{Files, Paths}

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.typesafe.config.ConfigFactory
import cromwell.CromwellSpec.IntegrationTest
import cromwell.engine.io.gcs.GoogleConfiguration
import cromwell.util.google.GoogleCredentialFactory.EnhancedCredentials
import org.scalatest.{Assertions, FlatSpec, Matchers}

import scala.util.Try

class GoogleCredentialFactorySpec extends FlatSpec with Matchers {

  behavior of "GoogleCredentialFactory"

  it should "refresh a token using user credentials" taggedAs IntegrationTest in {
    GoogleCredentialFactorySpec.assumeUserConfigExists()

    val credential = GoogleCredentialFactory.fromCromwellAuthScheme(GoogleCredentialFactorySpec.GoogleUserConfig)

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

  it should "refresh a token using a service account" taggedAs IntegrationTest in {
    GoogleCredentialFactorySpec.assumeAccountConfigExists()

    val credential = GoogleCredentialFactory.fromCromwellAuthScheme(GoogleCredentialFactorySpec.GoogleAccountConfig)

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

  it should "refresh a token using a refresh token" taggedAs IntegrationTest in {
    GoogleCredentialFactorySpec.assumeRefreshConfigExists()

    val refreshToken = GoogleCredentialFactory.fromCromwellAuthScheme(GoogleCredentialFactorySpec.GoogleRefreshConfig).freshCredential.get.getRefreshToken

    val firstUserCredentialsTry = GoogleCredentialFactory.fromUserAuthScheme(GoogleCredentialFactorySpec.GoogleRefreshConfig, refreshToken)

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
      .setTransport(GoogleCredentialFactory.httpTransport)
      .setJsonFactory(GoogleCredentialFactory.jsonFactory)
      .setClientSecrets("fakeId", "fakeSecret")
      .build()

    val exception = wrongCredentials.freshCredential.failed.get

    exception.getMessage should be("Unable to refresh token")
  }
}

object GoogleCredentialFactorySpec {
  val AccountConfigPath = Paths.get("cromwell-account.conf")
  val AccountConfigExists = Files.exists(AccountConfigPath)
  lazy val AccountConfig = ConfigFactory.parseFile(AccountConfigPath.toFile)
  lazy val GoogleAccountConfig: GoogleConfiguration = GoogleConfiguration.fromConfig(AccountConfig).get

  val UserConfigPath = Paths.get("cromwell-user.conf")
  val UserConfigExists = Files.exists(UserConfigPath)
  lazy val UserConfig = ConfigFactory.parseFile(UserConfigPath.toFile)
  lazy val GoogleUserConfig = GoogleConfiguration.fromConfig(UserConfig).get

  val RefreshConfigPath = Paths.get("cromwell-refresh.conf")
  val RefreshConfigExists = Files.exists(RefreshConfigPath)
  lazy val RefreshConfig = ConfigFactory.parseFile(RefreshConfigPath.toFile)
  lazy val GoogleRefreshConfig = GoogleConfiguration.fromConfig(RefreshConfig).get

  import Assertions._

  def assumeAccountConfigExists() = assume(AccountConfigExists, s"\nConfig not found $AccountConfigPath")

  def assumeUserConfigExists() = assume(UserConfigExists, s"\nConfig not found $UserConfigPath")

  def assumeRefreshConfigExists() = assume(RefreshConfigExists, s"\nConfig not found $RefreshConfigPath")
}
