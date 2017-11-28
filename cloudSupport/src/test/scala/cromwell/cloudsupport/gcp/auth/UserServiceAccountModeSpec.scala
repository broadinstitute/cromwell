package cromwell.cloudsupport.gcp.auth

import cromwell.cloudsupport.gcp.GoogleConfiguration
import org.scalatest.{FlatSpec, Matchers}

class UserServiceAccountModeSpec extends FlatSpec with Matchers {

  behavior of "UserServiceAccountMode"

  it should "generate a credential" in {
    val userServiceAccountMode = UserServiceAccountMode(
      "user-service-account",
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.userServiceAccountOptions
    val credentials = userServiceAccountMode.credential(workflowOptions)
    credentials.getAuthenticationType should be("OAuth2")
  }

  it should "validate with a user_service_account_json workflow option" in {
    val userServiceAccountMode = UserServiceAccountMode(
      "user-service-account",
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.userServiceAccountOptions
    userServiceAccountMode.validate(workflowOptions)
  }

  it should "fail validate without a user_service_account_json workflow option" in {
    val userServiceAccountMode = UserServiceAccountMode(
      "user-service-account",
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    val exception = intercept[OptionLookupException](userServiceAccountMode.validate(workflowOptions))
    exception.getMessage should be("user_service_account_json")
  }

  it should "requiresAuthFile" in {
    val userServiceAccountMode = UserServiceAccountMode(
      "user-service-account",
      GoogleConfiguration.GoogleScopes)
    userServiceAccountMode.requiresAuthFile should be(false)
  }

}
