package cromwell.cloudsupport.gcp.auth

import cromwell.core.TestKitSuite
import org.scalatest.{AsyncFlatSpecLike, Matchers}

class ApplicationDefaultModeSpec extends TestKitSuite("ApplicationDefaultModeSpec") with AsyncFlatSpecLike with Matchers {

  behavior of "ApplicationDefaultMode"

  it should "generate a credential" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()
    val applicationDefaultMode = new ApplicationDefaultMode("application-default")
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    applicationDefaultMode.credential(workflowOptions) map { credentials =>
      credentials.getAuthenticationType should be("OAuth2")
    }
  }

  it should "validate" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()
    val applicationDefaultMode = new ApplicationDefaultMode("application-default")
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    applicationDefaultMode.validate(workflowOptions)
    succeed
  }

  it should "requiresAuthFile" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()
    val applicationDefaultMode = new ApplicationDefaultMode("application-default")
    applicationDefaultMode.requiresAuthFile should be(false)
    succeed
  }

}
