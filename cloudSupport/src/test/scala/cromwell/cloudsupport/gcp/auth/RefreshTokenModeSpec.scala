package cromwell.cloudsupport.gcp.auth

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class RefreshTokenModeSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "RefreshTokenMode"

  it should "fail to generate a bad credential" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret"
    )
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    val exception = intercept[RuntimeException](refreshTokenMode.credentials(workflowOptions))
    exception.getMessage should startWith("Google credentials are invalid: ")
  }

  it should "fail to generate a credential that cannot be validated" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret"
    )
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    refreshTokenMode.credentialsValidation = _ => throw new IllegalArgumentException("no worries! this is expected")
    val exception = intercept[RuntimeException](refreshTokenMode.credentials(workflowOptions))
    exception.getMessage should startWith("Google credentials are invalid: ")
  }

  it should "generate a non-validated credential" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret"
    )
    refreshTokenMode.credentialsValidation = GoogleAuthMode.NoCredentialsValidation
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    val credentials = refreshTokenMode.credentials(workflowOptions)
    credentials.getAuthenticationType should be("OAuth2")
  }

  it should "fail to generate without a refresh_token workflow option" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret"
    )
    val exception = intercept[OptionLookupException](refreshTokenMode.credentials())
    exception.getMessage should be("refresh_token")
  }

  it should "requiresAuthFile" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret"
    )
    refreshTokenMode.requiresAuthFile should be(true)
  }

}
