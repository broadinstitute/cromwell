package cromwell.cloudsupport.gcp.auth

import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.{CromwellFatalException, TestKitSuite}
import org.scalatest.{AsyncFlatSpecLike, Matchers}

class RefreshTokenModeSpec extends TestKitSuite("RefreshTokenModeSpec") with AsyncFlatSpecLike with Matchers {

  behavior of "RefreshTokenMode"

  it should "fail to generate a bad credential" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret",
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    recoverToExceptionIf[CromwellFatalException](refreshTokenMode.credential(workflowOptions)) map { exception =>
      exception.getCause.getMessage should startWith("Google credentials are invalid: ")
    }
  }

  it should "fail to generate a credential that cannot be validated" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret",
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    refreshTokenMode.credentialValidation = _ => throw new IllegalArgumentException("no worries! this is expected")
    recoverToExceptionIf[CromwellFatalException](refreshTokenMode.credential(workflowOptions)) map { exception =>
      exception.getCause.getMessage should startWith("Google credentials are invalid: ")
    }
  }

  it should "generate a non-validated credential" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret",
      GoogleConfiguration.GoogleScopes)
    refreshTokenMode.credentialValidation = _ => ()
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    refreshTokenMode.credential(workflowOptions) map { credentials =>
      credentials.getAuthenticationType should be("OAuth2")
    }
  }

  it should "validate with a refresh_token workflow option" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret",
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.refreshTokenOptions
    refreshTokenMode.validate(workflowOptions)
    succeed
  }

  it should "fail validate without a refresh_token workflow option" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret",
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    val exception = intercept[IllegalArgumentException](refreshTokenMode.validate(workflowOptions))
    exception.getMessage should be("Missing parameters in workflow options: refresh_token")
    succeed
  }

  it should "requiresAuthFile" in {
    val refreshTokenMode = RefreshTokenMode(
      "user-via-refresh",
      "secret_id",
      "secret_secret",
      GoogleConfiguration.GoogleScopes)
    refreshTokenMode.requiresAuthFile should be(true)
    succeed
  }

}
