package cromwell.cloudsupport.gcp.auth

import java.io.FileNotFoundException

import better.files.File
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.TestKitSuite
import org.scalatest.{AsyncFlatSpecLike, Matchers}

class ServiceAccountModeSpec extends TestKitSuite("ServiceAccountModeSpec") with AsyncFlatSpecLike with Matchers {

  behavior of "ServiceAccountMode"

  it should "fail to generate a bad credential from json" in {
    val jsonMockFile = File
      .newTemporaryFile("service-account.", ".json")
      .write(GoogleAuthModeSpec.serviceAccountJsonContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.JsonFileFormat(jsonMockFile.pathAsString),
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    val exception = intercept[RuntimeException](serviceAccountMode.credential(workflowOptions))
    exception.getMessage should startWith("Google credentials are invalid: ")
    jsonMockFile.delete(true)
    succeed
  }

  it should "fail to generate a bad credential from a pem" in {
    val pemMockFile = File
      .newTemporaryFile("service-account.", ".pem")
      .write(GoogleAuthModeSpec.serviceAccountPemContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.PemFileFormat("the_account_id", pemMockFile.pathAsString),
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    val exception = intercept[RuntimeException](serviceAccountMode.credential(workflowOptions))
    exception.getMessage should startWith("Google credentials are invalid: ")
    pemMockFile.delete(true)
    succeed
  }

  it should "fail to generate a bad credential from a missing json" in {
    val jsonMockFile = File.newTemporaryFile("service-account.", ".json").delete()
    val exception = intercept[FileNotFoundException] {
      ServiceAccountMode(
        "service-account",
        ServiceAccountMode.JsonFileFormat(jsonMockFile.pathAsString),
        GoogleConfiguration.GoogleScopes)
    }
    exception.getMessage should fullyMatch regex "File .*/service-account..*.json does not exist or is not readable"
    succeed
  }

  it should "fail to generate a bad credential from a missing pem" in {
    val pemMockFile = File.newTemporaryFile("service-account.", ".pem").delete()
    val exception = intercept[FileNotFoundException] {
      ServiceAccountMode(
        "service-account",
        ServiceAccountMode.PemFileFormat("the_account_id", pemMockFile.pathAsString),
        GoogleConfiguration.GoogleScopes)
    }
    exception.getMessage should fullyMatch regex "File .*/service-account..*.pem does not exist or is not readable"
    succeed
  }

  it should "generate a non-validated credential from json" in {
    val jsonMockFile = File
      .newTemporaryFile("service-account.", ".json")
      .write(GoogleAuthModeSpec.serviceAccountJsonContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.JsonFileFormat(jsonMockFile.pathAsString),
      GoogleConfiguration.GoogleScopes)
    serviceAccountMode.credentialValidation = _ => ()
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    serviceAccountMode.credential(workflowOptions) map { credentials =>
      credentials.getAuthenticationType should be("OAuth2")
      jsonMockFile.delete(true)
      succeed
    }
  }

  it should "generate a non-validated credential from a pem" in {
    val pemMockFile = File
      .newTemporaryFile("service-account.", ".pem")
      .write(GoogleAuthModeSpec.serviceAccountPemContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.PemFileFormat("the_account_id", pemMockFile.pathAsString),
      GoogleConfiguration.GoogleScopes)
    serviceAccountMode.credentialValidation = _ => ()
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    serviceAccountMode.credential(workflowOptions) map { credentials =>
      credentials.getAuthenticationType should be("OAuth2")
      pemMockFile.delete(true)
      succeed
    }
  }

  it should "pass validate with a refresh_token workflow option from json" in {
    val jsonMockFile = File
      .newTemporaryFile("service-account.", ".json")
      .write(GoogleAuthModeSpec.serviceAccountJsonContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.JsonFileFormat(jsonMockFile.pathAsString),
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    serviceAccountMode.validate(workflowOptions)
    jsonMockFile.delete(true)
    succeed
  }

  it should "pass validate with a refresh_token workflow option from a pem" in {
    val pemMockFile = File
      .newTemporaryFile("service-account.", ".pem")
      .write(GoogleAuthModeSpec.serviceAccountPemContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.PemFileFormat("the_account_id", pemMockFile.pathAsString),
      GoogleConfiguration.GoogleScopes)
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    serviceAccountMode.validate(workflowOptions)
    pemMockFile.delete(true)
    succeed
  }

  it should "requiresAuthFile from json" in {
    val jsonMockFile = File
      .newTemporaryFile("service-account.", ".json")
      .write(GoogleAuthModeSpec.serviceAccountJsonContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.JsonFileFormat(jsonMockFile.pathAsString),
      GoogleConfiguration.GoogleScopes)
    serviceAccountMode.requiresAuthFile should be(false)
    jsonMockFile.delete(true)
    succeed
  }

  it should "requiresAuthFile from a pem" in {
    val pemMockFile = File
      .newTemporaryFile("service-account.", ".pem")
      .write(GoogleAuthModeSpec.serviceAccountPemContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.PemFileFormat("the_account_id", pemMockFile.pathAsString),
      GoogleConfiguration.GoogleScopes)
    serviceAccountMode.requiresAuthFile should be(false)
    pemMockFile.delete(true)
    succeed
  }

}
