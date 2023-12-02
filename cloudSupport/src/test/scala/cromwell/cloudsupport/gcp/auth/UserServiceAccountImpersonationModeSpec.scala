package cromwell.cloudsupport.gcp.auth

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class UserServiceAccountImpersonationModeSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "UserServiceAccountImpersonationMode"

  it should "generate a non-validated credential" in {
    val impersonationMode = UserServiceAccountImpersonationMode("user-service-account-impersonation")
    val workflowOptions = GoogleAuthModeSpec.userServiceAccountImpersonationOptions
    impersonationMode.credentialsValidation = GoogleAuthMode.NoCredentialsValidation
    val credentials = impersonationMode.credentials(workflowOptions)
    credentials.getAuthenticationType should be("OAuth2")
  }

  it should "fail to generate credentials without a user_service_account_email workflow option" in {
    val impersonationMode = UserServiceAccountImpersonationMode("user-service-account-impersonation")
    val exception = intercept[OptionLookupException](impersonationMode.credentials())
    exception.getMessage should be("user_service_account_email")
  }
}
