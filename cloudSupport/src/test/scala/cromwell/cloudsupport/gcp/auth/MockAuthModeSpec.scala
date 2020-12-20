package cromwell.cloudsupport.gcp.auth

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class MockAuthModeSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "MockAuthMode"

  it should "generate a credential" in {
    val mockAuthMode = MockAuthMode("no_auth")
    val credentials = mockAuthMode.credentials()
    credentials.getAuthenticationType should be("OAuth2")
  }

  it should "requiresAuthFile" in {
    val mockAuthMode = MockAuthMode("no_auth")
    mockAuthMode.requiresAuthFile should be(false)
    succeed
  }

}
