package cromwell.util.google

import java.nio.file.{Files, Paths}

import com.google.api.client.googleapis.auth.oauth2.GoogleCredential
import com.typesafe.config.ConfigFactory
import org.scalatest.{Assertions, FlatSpec, Matchers}

import scala.util.Try

class GoogleCredentialFactorySpec extends FlatSpec with Matchers {

  behavior of "GoogleCredentialFactory"

  it should "refresh a token using the default application credential" in {
    GoogleCredentialFactorySpec.assumeApplicationDefaultCredentialExists()

    val credentialFactory = new ApplicationDefaultCredentialFactory()
    testFactory(credentialFactory)
  }

  it should "refresh a token using user credentials" in {
    GoogleCredentialFactorySpec.assumeUserConfigExists()

    val credentialFactory = GoogleCredentialFactory.fromAuthScheme(GoogleCredentialFactorySpec.UserConfig)
    testFactory(credentialFactory)
  }

  it should "refresh a token using a service account" in {
    GoogleCredentialFactorySpec.assumeAccountConfigExists()

    val credentialFactory = GoogleCredentialFactory.fromAuthScheme(GoogleCredentialFactorySpec.AccountConfig)
    testFactory(credentialFactory)
  }

  it should "refresh a token using a refresh token" in {
    GoogleCredentialFactorySpec.assumeRefreshConfigExists()

    val refreshTokenFactory = GoogleCredentialFactory.fromAuthScheme(GoogleCredentialFactorySpec.RefreshConfig)
    val refreshToken = refreshTokenFactory.freshCredential.get.getRefreshToken

    val credentialFactory = new RefreshTokenCredentialFactory(GoogleCredentialFactorySpec.RefreshConfig, refreshToken)
    testFactory(credentialFactory)
  }

  it should "not refresh an empty token" in {
    val factory = new GoogleCredentialFactory {
      override protected def initCredential() = new GoogleCredential()
    }

    val exception = factory.freshCredential.failed.get
    exception.getMessage should be("Unable to refresh token")
  }

  private[this] def testFactory(credentialFactory: GoogleCredentialFactory) {
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
}

object GoogleCredentialFactorySpec {
  val ApplicationDefaultCredentialTry = Try(new ApplicationDefaultCredentialFactory().rawCredential)
  lazy val ApplicationDefaultCredentialExists = ApplicationDefaultCredentialTry.isSuccess

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

  def assumeApplicationDefaultCredentialExists() = assume(ApplicationDefaultCredentialExists,
    // Use fold to defer execution of the failed message.
    ApplicationDefaultCredentialTry.toOption.fold(
      s"""
         |Error getting default application credential. Please see this link for more info:
         |https://developers.google.com/accounts/docs/application-default-credentials
         |
         |${ApplicationDefaultCredentialTry.failed.get.getMessage}""".stripMargin)(_ => "aok"))

  def assumeAccountConfigExists() = assume(AccountConfigExists, s"\nConfig not found $AccountConfigPath")

  def assumeUserConfigExists() = assume(UserConfigExists, s"\nConfig not found $UserConfigPath")

  def assumeRefreshConfigExists() = assume(RefreshConfigExists, s"\nConfig not found $RefreshConfigPath")
}
