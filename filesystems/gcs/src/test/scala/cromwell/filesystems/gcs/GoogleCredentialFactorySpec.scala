package cromwell.filesystems.gcs

import java.nio.file.{Files, Paths}

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.typesafe.config.ConfigFactory
import cromwell.filesystems.gcs.GcsFileSystemSpec.IntegrationTest
import cromwell.filesystems.gcs.GoogleCredentialFactory.EnhancedCredentials
import org.scalatest.{Assertions, FlatSpec, Matchers}

import scala.util.Try

class GoogleCredentialFactorySpec extends FlatSpec with Matchers {

  behavior of "GoogleCredentialFactory"

  it should "refresh a token using user credentials" taggedAs IntegrationTest in {
    GoogleCredentialFactorySpec.assumeUserConfigExists()

    val credentialFactory = GoogleCredentialFactory(GoogleCredentialFactorySpec.GoogleUserConfig.authMode, GcsScopes)

    val firstCredentialTry = credentialFactory.freshCredential
    assert(firstCredentialTry.isSuccess)
    val firstCredential = firstCredentialTry.get
    firstCredential.getAccessToken shouldNot be(empty)

    firstCredential.setExpiresInSeconds(59L)

    val secondCredentialTry = credentialFactory.freshCredential
    assert(secondCredentialTry.isSuccess)

    val secondCredential = secondCredentialTry.get
    secondCredential.getAccessToken shouldNot be(empty)
    secondCredential.getExpiresInSeconds shouldNot be(null)
    secondCredential.getExpiresInSeconds.longValue should be > 60L
  }

  it should "refresh a token using a service account" taggedAs IntegrationTest in {
    GoogleCredentialFactorySpec.assumeAccountConfigExists()

    val credentialFactory = GoogleCredentialFactory(GoogleCredentialFactorySpec.GoogleAccountConfig.authMode, GcsScopes)

    val firstCredentialTry = credentialFactory.freshCredential
    assert(firstCredentialTry.isSuccess)
    val firstCredential = firstCredentialTry.get
    firstCredential.getAccessToken shouldNot be(empty)

    firstCredential.setExpiresInSeconds(59L)

    val secondCredentialTry = credentialFactory.freshCredential
    assert(secondCredentialTry.isSuccess)

    val secondCredential = secondCredentialTry.get
    secondCredential.getAccessToken shouldNot be(empty)
    secondCredential.getExpiresInSeconds shouldNot be(null)
    secondCredential.getExpiresInSeconds.longValue should be > 60L
  }

  it should "refresh a token using a refresh token" taggedAs IntegrationTest in {
    GoogleCredentialFactorySpec.assumeRefreshConfigExists()

    val refreshToken = {
      val googleConf = GoogleConfigurationAdapter.build(GoogleCredentialFactorySpec.RefreshConfig).cromwellConf
      GoogleCredentialFactory(googleConf.authMode, GcsScopes).freshCredential.get.getRefreshToken
    }

    val credentialFactoryBuilder =  Try {
      GoogleCredentialFactory(GoogleCredentialFactorySpec.GoogleAccountConfig.authMode, GcsScopes, Option(refreshToken))
    }

    assert(credentialFactoryBuilder.isSuccess)
    val credentialFactory = credentialFactoryBuilder.get

    val firstCredentialTry = credentialFactory.freshCredential
    assert(firstCredentialTry.isSuccess)
    val firstCredential = firstCredentialTry.get
    firstCredential.getAccessToken shouldNot be(empty)

    firstCredential.setExpiresInSeconds(59L)

    val secondCredentialTry = credentialFactory.freshCredential
    assert(secondCredentialTry.isSuccess)

    val secondCredential = secondCredentialTry.get
    secondCredential.getAccessToken shouldNot be(empty)
    secondCredential.getExpiresInSeconds shouldNot be(null)
    secondCredential.getExpiresInSeconds.longValue should be > 60L
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
  lazy val GoogleAccountConfig = GoogleConfigurationAdapter.build(AccountConfig).cromwellConf

  val UserConfigPath = Paths.get("cromwell-user.conf")
  val UserConfigExists = Files.exists(UserConfigPath)
  lazy val UserConfig = ConfigFactory.parseFile(UserConfigPath.toFile)
  lazy val GoogleUserConfig = GoogleConfigurationAdapter.build(UserConfig).cromwellConf

  val RefreshConfigPath = Paths.get("cromwell-refresh.conf")
  val RefreshConfigExists = Files.exists(RefreshConfigPath)
  lazy val RefreshConfig = ConfigFactory.parseFile(RefreshConfigPath.toFile)
  lazy val GoogleRefreshConfig = GoogleConfigurationAdapter.build(RefreshConfig).userConf.get

  import Assertions._

  def assumeAccountConfigExists() = assume(AccountConfigExists, s"\nConfig not found $AccountConfigPath")

  def assumeUserConfigExists() = assume(UserConfigExists, s"\nConfig not found $UserConfigPath")

  def assumeRefreshConfigExists() = assume(RefreshConfigExists, s"\nConfig not found $RefreshConfigPath")
}
