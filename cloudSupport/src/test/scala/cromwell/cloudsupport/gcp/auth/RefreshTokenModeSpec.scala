package cromwell.cloudsupport.gcp.auth

import org.scalatest.{FlatSpec, Matchers}

class RefreshTokenModeSpec extends FlatSpec with Matchers {

  behavior of "RefreshTokenMode"

  it should "fail to generate a bad credential" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret")
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    val exception = intercept[RuntimeException](refreshTokenMode.credential(workflowOptions))
    exception.getMessage should startWith("Google credentials are invalid: ")
  }

  it should "fail to generate a credential that cannot be validated" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret")
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    refreshTokenMode.credentialValidation = _ => throw new IllegalArgumentException("no worries! this is expected")
    val exception = intercept[RuntimeException](refreshTokenMode.credential(workflowOptions))
    exception.getMessage should startWith("Google credentials are invalid: ")
  }

  it should "generate a non-validated credential" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret")
    refreshTokenMode.credentialValidation = _ => ()
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    val credentials = refreshTokenMode.credential(workflowOptions)
    credentials.getAuthenticationType should be("OAuth2")
  }

  it should "validate with a refresh_token workflow option" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret")
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    refreshTokenMode.validate(workflowOptions)
  }

  it should "fail validate without a refresh_token workflow option" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret")
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    val exception = intercept[OptionLookupException](refreshTokenMode.validate(workflowOptions))
    exception.getMessage should be("refresh_token")
  }

  it should "requiresAuthFile" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret")
    refreshTokenMode.requiresAuthFile should be(true)
  }

}
