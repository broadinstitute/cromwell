package cromwell.util.google

import java.nio.file.{Files, Paths}

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.typesafe.config.ConfigFactory
import cromwell.engine.io.gcs.GoogleConfiguration
import cromwell.util.google.GoogleCredentialFactory.EnhancedCredentials
import org.scalatest.{Assertions, FlatSpec, Matchers}

class GoogleCredentialFactorySpec extends FlatSpec with Matchers {

  behavior of "GoogleCredentialFactory"

  it should "refresh a token using user credentials" in {
    GoogleCredentialFactorySpec.assumeUserConfigExists()

    val credentialFactory = {
      new GoogleCredentialFactory() {
        override val GoogleConf = GoogleCredentialFactorySpec.GoogleUserConfig
      }.fromCromwellAuthScheme
    }

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

  it should "refresh a token using a service account" in {
    GoogleCredentialFactorySpec.assumeAccountConfigExists()

    val credentialFactory = {
      new GoogleCredentialFactory() {
        override val GoogleConf = GoogleCredentialFactorySpec.GoogleAccountConfig
      }.fromCromwellAuthScheme
    }

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

  it should "refresh a token using a refresh token" in {
    GoogleCredentialFactorySpec.assumeRefreshConfigExists()

    val refreshToken = {
      new GoogleCredentialFactory() {
        override val GoogleConf = GoogleCredentialFactorySpec.GoogleRefreshConfig
      }.fromCromwellAuthScheme.freshCredential.get.getRefreshToken
    }

    val credentialFactoryBuilder = new GoogleCredentialFactory() {
      override val GoogleConf = GoogleCredentialFactorySpec.GoogleRefreshConfig
    }.fromUserAuthScheme(refreshToken)

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
  lazy val GoogleAccountConfig = GoogleConfiguration.build(AccountConfig)

  val UserConfigPath = Paths.get("cromwell-user.conf")
  val UserConfigExists = Files.exists(UserConfigPath)
  lazy val UserConfig = ConfigFactory.parseFile(UserConfigPath.toFile)
  lazy val GoogleUserConfig = GoogleConfiguration.build(UserConfig)

  val RefreshConfigPath = Paths.get("cromwell-refresh.conf")
  val RefreshConfigExists = Files.exists(RefreshConfigPath)
  lazy val RefreshConfig = ConfigFactory.parseFile(RefreshConfigPath.toFile)
  lazy val GoogleRefreshConfig = GoogleConfiguration.build(RefreshConfig)

  import Assertions._

  def assumeAccountConfigExists() = assume(AccountConfigExists, s"\nConfig not found $AccountConfigPath")

  def assumeUserConfigExists() = assume(UserConfigExists, s"\nConfig not found $UserConfigPath")

  def assumeRefreshConfigExists() = assume(RefreshConfigExists, s"\nConfig not found $RefreshConfigPath")
}
