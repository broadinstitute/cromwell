package cromwell.cloudsupport.gcp.auth

import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.TestKitSuite
import org.scalatest.{AsyncFlatSpecLike, Matchers}

class UserServiceAccountModeSpec extends TestKitSuite("UserServiceAccountModeSpec") with AsyncFlatSpecLike with Matchers {

  behavior of "UserServiceAccountMode"

  it should "generate a credential" in {
    val userServiceAccountMode = UserServiceAccountMode(
      "user-service-account",
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.userServiceAccountOptions
    userServiceAccountMode.credential(workflowOptions) map { credentials =>
      credentials.getAuthenticationType should be("OAuth2")
    }
  }

  it should "validate with a user_service_account_json workflow option" in {
    val userServiceAccountMode = UserServiceAccountMode(
      "user-service-account",
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.userServiceAccountOptions
    userServiceAccountMode.validate(workflowOptions)
    succeed
  }

  it should "fail validate without a user_service_account_json workflow option" in {
    val userServiceAccountMode = UserServiceAccountMode(
      "user-service-account",
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    val exception = intercept[IllegalArgumentException](userServiceAccountMode.validate(workflowOptions))
    exception.getMessage should be("Missing parameters in workflow options: user_service_account_json")
    succeed
  }

  it should "requiresAuthFile" in {
    val userServiceAccountMode = UserServiceAccountMode(
      "user-service-account",
      GoogleConfiguration.GoogleScopes)
    userServiceAccountMode.requiresAuthFile should be(false)
    succeed
  }

}
