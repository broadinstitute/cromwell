package cromwell.cloudsupport.gcp.auth

import org.scalatest.{FlatSpec, Matchers}

class MockAuthModeSpec extends FlatSpec with Matchers {

  behavior of "MockAuthMode"

  it should "generate a credential" in {
    val mockAuthMode = MockAuthMode
    val credentials = mockAuthMode.credentials()
    credentials.getAuthenticationType should be("OAuth2")
  }

  it should "requiresAuthFile" in {
    val mockAuthMode = MockAuthMode
    mockAuthMode.requiresAuthFile should be(false)
    succeed
  }

}
