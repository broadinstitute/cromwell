package cromwell.filesystems.gcs

import java.nio.file.{Files, Paths}

import com.google.api.client.auth.oauth2.Credential
import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.jackson2.JacksonFactory
import com.typesafe.config.ConfigFactory
import cromwell.filesystems.gcs.GcsFileSystemSpec.IntegrationTest
import cromwell.filesystems.gcs.GoogleAuthMode.EnhancedCredentials
import org.scalatest.{Assertions, FlatSpec, Matchers}

import scala.util.Try


class GoogleCredentialFactorySpec extends FlatSpec with Matchers {
  import GoogleCredentialFactorySpec._

  val options = GoogleOptionsMap(Map.empty)

  behavior of "GoogleCredentialFactory"

  it should "refresh a token using user credentials" taggedAs IntegrationTest in {
    GoogleCredentialFactorySpec.assumeUserConfigExists()

    def user(k: String): String = UserConfig.getString(k)

    val credential = UserMode(
      name = "user",
      user = user("user"),
      secretsFile = user("secrets-file"),
      datastoreDir = user("data-store-dir")).credential(options)

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

    val credential = ApplicationDefaultMode(name = "default").credential(options)

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

    def refresh(k: String): String = AccountConfig.getString(k)
    val opts = GoogleOptionsMap(Map("refresh_token" -> refresh("refresh_token")))

    val credential = RefreshTokenMode(name = "refresh", refresh("client-id"), refresh("client-secret")).credential(opts)

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
}

object GoogleCredentialFactorySpec {
  val AccountConfigPath = Paths.get("cromwell-account.conf")
  val AccountConfigExists = Files.exists(AccountConfigPath)
  lazy val AccountConfig = ConfigFactory.parseFile(AccountConfigPath.toFile)

  val UserConfigPath = Paths.get("cromwell-user.conf")
  val UserConfigExists = Files.exists(UserConfigPath)
  lazy val UserConfig = ConfigFactory.parseFile(UserConfigPath.toFile)

  val RefreshConfigPath = Paths.get("cromwell-refresh.conf")
  val RefreshConfigExists = Files.exists(RefreshConfigPath)
  lazy val RefreshConfig = ConfigFactory.parseFile(RefreshConfigPath.toFile)

  import Assertions._

  def assumeAccountConfigExists() = assume(AccountConfigExists, s"\nConfig not found $AccountConfigPath")

  def assumeUserConfigExists() = assume(UserConfigExists, s"\nConfig not found $UserConfigPath")

  def assumeRefreshConfigExists() = assume(RefreshConfigExists, s"\nConfig not found $RefreshConfigPath")

  case class GoogleOptionsMap(map: Map[String, String]) extends GoogleAuthMode.GoogleAuthOptions {
    override def get(key: String): Try[String] = Try { map.get(key).get }
  }
}
