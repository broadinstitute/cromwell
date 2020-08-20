package cromwell.cloudsupport.gcp.auth

import java.io.FileNotFoundException

import better.files.File
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class ServiceAccountModeSpec extends AnyFlatSpec with Matchers {

  behavior of "ServiceAccountMode"

  it should "fail to generate a bad credential from json" in {
    val jsonMockFile = File
      .newTemporaryFile("service-account.", ".json")
      .write(GoogleAuthModeSpec.serviceAccountJsonContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.JsonFileFormat(jsonMockFile.pathAsString)
    )
    val exception = intercept[RuntimeException](serviceAccountMode.credentials())
    exception.getMessage should startWith("Google credentials are invalid: ")
    jsonMockFile.delete(true)
  }

  it should "fail to generate a bad credential from a pem" in {
    val pemMockFile = File
      .newTemporaryFile("service-account.", ".pem")
      .write(GoogleAuthModeSpec.serviceAccountPemContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.PemFileFormat("the_account_id", pemMockFile.pathAsString),
    )
    val exception = intercept[RuntimeException](serviceAccountMode.credentials())
    exception.getMessage should startWith("Google credentials are invalid: ")
    pemMockFile.delete(true)
  }

  it should "fail to generate a bad credential from a missing json" in {
    val jsonMockFile = File.newTemporaryFile("service-account.", ".json").delete()
    val exception = intercept[FileNotFoundException] {
      ServiceAccountMode(
        "service-account",
        ServiceAccountMode.JsonFileFormat(jsonMockFile.pathAsString)
      )
    }
    exception.getMessage should fullyMatch regex "File .*/service-account..*.json does not exist or is not readable"
  }

  it should "fail to generate a bad credential from a missing pem" in {
    val pemMockFile = File.newTemporaryFile("service-account.", ".pem").delete()
    val exception = intercept[FileNotFoundException] {
      ServiceAccountMode(
        "service-account",
        ServiceAccountMode.PemFileFormat("the_account_id", pemMockFile.pathAsString),
      )
    }
    exception.getMessage should fullyMatch regex "File .*/service-account..*.pem does not exist or is not readable"
  }

  it should "generate a non-validated credential from json" in {
    val jsonMockFile = File
      .newTemporaryFile("service-account.", ".json")
      .write(GoogleAuthModeSpec.serviceAccountJsonContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.JsonFileFormat(jsonMockFile.pathAsString)
    )
    serviceAccountMode.credentialsValidation = GoogleAuthMode.NoCredentialsValidation
    val credentials = serviceAccountMode.credentials()
    credentials.getAuthenticationType should be("OAuth2")
    jsonMockFile.delete(true)
  }

  it should "generate a non-validated credential from a pem" in {
    val pemMockFile = File
      .newTemporaryFile("service-account.", ".pem")
      .write(GoogleAuthModeSpec.serviceAccountPemContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.PemFileFormat("the_account_id", pemMockFile.pathAsString),
    )
    serviceAccountMode.credentialsValidation = GoogleAuthMode.NoCredentialsValidation
    val credentials = serviceAccountMode.credentials()
    credentials.getAuthenticationType should be("OAuth2")
    pemMockFile.delete(true)
  }

  it should "requiresAuthFile from json" in {
    val jsonMockFile = File
      .newTemporaryFile("service-account.", ".json")
      .write(GoogleAuthModeSpec.serviceAccountJsonContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.JsonFileFormat(jsonMockFile.pathAsString)
    )
    serviceAccountMode.requiresAuthFile should be(false)
    jsonMockFile.delete(true)
  }

  it should "requiresAuthFile from a pem" in {
    val pemMockFile = File
      .newTemporaryFile("service-account.", ".pem")
      .write(GoogleAuthModeSpec.serviceAccountPemContents)
    val serviceAccountMode = ServiceAccountMode(
      "service-account",
      ServiceAccountMode.PemFileFormat("the_account_id", pemMockFile.pathAsString)
    )
    serviceAccountMode.requiresAuthFile should be(false)
    pemMockFile.delete(true)
  }

}
