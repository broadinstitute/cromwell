package cromwell.cloudsupport.gcp.auth

import org.scalatest.{FlatSpec, Matchers}

class ApplicationDefaultModeSpec extends FlatSpec with Matchers {

  behavior of "ApplicationDefaultMode"

  it should "generate a credential" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()
    val applicationDefaultMode = new ApplicationDefaultMode("application-default")
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    val credentials = applicationDefaultMode.credential(workflowOptions)
    credentials.getAuthenticationType should be("OAuth2")
  }

  it should "validate" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()
    val applicationDefaultMode = new ApplicationDefaultMode("application-default")
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    applicationDefaultMode.validate(workflowOptions)
  }

  it should "requiresAuthFile" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()
    val applicationDefaultMode = new ApplicationDefaultMode("application-default")
    applicationDefaultMode.requiresAuthFile should be(false)
  }

}
