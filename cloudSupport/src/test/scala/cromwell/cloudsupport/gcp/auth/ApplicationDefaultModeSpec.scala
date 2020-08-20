package cromwell.cloudsupport.gcp.auth

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class ApplicationDefaultModeSpec extends AnyFlatSpec with Matchers {

  behavior of "ApplicationDefaultMode"

  it should "generate a credential" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()
    val applicationDefaultMode = new ApplicationDefaultMode("application-default")
    val credentials = applicationDefaultMode.credentials()
    credentials.getAuthenticationType should be("OAuth2")
  }

  it should "requiresAuthFile" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()
    val applicationDefaultMode = new ApplicationDefaultMode("application-default")
    applicationDefaultMode.requiresAuthFile should be(false)
  }

}
