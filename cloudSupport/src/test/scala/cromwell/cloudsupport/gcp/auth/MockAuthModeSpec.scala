package cromwell.cloudsupport.gcp.auth

import cromwell.core.TestKitSuite
import org.scalatest.{AsyncFlatSpecLike, Matchers}

class MockAuthModeSpec extends TestKitSuite("MockAuthModeSpec") with AsyncFlatSpecLike with Matchers {

  behavior of "MockAuthMode"

  it should "generate a credential" in {
    val mockAuthMode = MockAuthMode
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    mockAuthMode.credential(workflowOptions) map { credentials =>
      credentials.getAuthenticationType should be("OAuth2")
    }
  }

  it should "validate" in {
    val mockAuthMode = MockAuthMode
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    mockAuthMode.validate(workflowOptions)
    succeed
  }

  it should "requiresAuthFile" in {
    val mockAuthMode = MockAuthMode
    mockAuthMode.requiresAuthFile should be(false)
    succeed
  }

}
