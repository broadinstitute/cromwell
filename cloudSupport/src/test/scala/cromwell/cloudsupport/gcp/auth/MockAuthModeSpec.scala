package cromwell.cloudsupport.gcp.auth

import org.scalatest.{FlatSpec, Matchers}

class MockAuthModeSpec extends FlatSpec with Matchers {

  behavior of "MockAuthMode"

  it should "generate a credential" in {
    val mockAuthMode = MockAuthMode
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    val credentials = mockAuthMode.credential(workflowOptions)
    credentials.getAuthenticationType should be("OAuth2")
  }

  it should "validate" in {
    val mockAuthMode = MockAuthMode
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    mockAuthMode.validate(workflowOptions)
  }

  it should "requiresAuthFile" in {
    val mockAuthMode = MockAuthMode
    mockAuthMode.requiresAuthFile should be(false)
    succeed
  }

}
