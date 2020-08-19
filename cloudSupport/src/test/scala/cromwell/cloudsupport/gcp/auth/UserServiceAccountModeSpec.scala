package cromwell.cloudsupport.gcp.auth

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class UserServiceAccountModeSpec extends AnyFlatSpec with Matchers {

  behavior of "UserServiceAccountMode"

  it should "generate a non-validated credential" in {
    val userServiceAccountMode = UserServiceAccountMode("user-service-account")
    val workflowOptions = GoogleAuthModeSpec.userServiceAccountOptions
    userServiceAccountMode.credentialsValidation = GoogleAuthMode.NoCredentialsValidation
    val credentials = userServiceAccountMode.credentials(workflowOptions)
    credentials.getAuthenticationType should be("OAuth2")
  }

  it should "fail to generate credentials without a user_service_account_json workflow option" in {
    val userServiceAccountMode = UserServiceAccountMode("user-service-account")
    val exception = intercept[OptionLookupException](userServiceAccountMode.credentials())
    exception.getMessage should be("user_service_account_json")
  }

  it should "requiresAuthFile" in {
    val userServiceAccountMode = UserServiceAccountMode("user-service-account")
    userServiceAccountMode.requiresAuthFile should be(false)
  }

}
