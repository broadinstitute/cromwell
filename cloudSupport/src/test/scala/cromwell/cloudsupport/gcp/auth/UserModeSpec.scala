package cromwell.cloudsupport.gcp.auth

import java.io.FileNotFoundException

import better.files.File
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.TestKitSuite
import org.scalatest.{AsyncFlatSpecLike, Matchers}

class UserModeSpec extends TestKitSuite("UserModeSpec") with AsyncFlatSpecLike with Matchers {

  behavior of "UserMode"

  it should "fail to generate a bad credential" in {
    val secretsMockFile = File.newTemporaryFile("secrets.", ".json").write(GoogleAuthModeSpec.userCredentialsContents)
    val dataStoreMockDir = File.newTemporaryDirectory("dataStore.")
    val userMode = UserMode(
      "user",
      "alice",
      secretsMockFile.pathAsString,
      dataStoreMockDir.pathAsString,
      GoogleConfiguration.GoogleScopes
    )
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    val exception = intercept[RuntimeException](userMode.credential(workflowOptions))
    exception.getMessage should startWith("Google credentials are invalid: ")
    secretsMockFile.delete(true)
    dataStoreMockDir.delete(true)
    succeed
  }

  it should "fail to generate a bad credential from a secrets json" in {
    val secretsMockFile = File.newTemporaryFile("secrets.", ".json").delete()
    val dataStoreMockDir = File.newTemporaryDirectory("dataStore.")
    val userMode = UserMode(
      "user",
      "alice",
      secretsMockFile.pathAsString,
      dataStoreMockDir.pathAsString,
      GoogleConfiguration.GoogleScopes
    )
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    val exception = intercept[FileNotFoundException](userMode.credential(workflowOptions))
    exception.getMessage should fullyMatch regex "File .*/secrets..*.json does not exist or is not readable"
    dataStoreMockDir.delete(true)
    succeed
  }

  it should "generate a non-validated credential" in {
    val secretsMockFile = File.newTemporaryFile("secrets.", ".json").write(GoogleAuthModeSpec.userCredentialsContents)
    val dataStoreMockDir = File.newTemporaryDirectory("dataStore.")
    val userMode = UserMode(
      "user",
      "alice",
      secretsMockFile.pathAsString,
      dataStoreMockDir.pathAsString,
      GoogleConfiguration.GoogleScopes
    )
    userMode.credentialValidation = _ => ()
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    userMode.credential(workflowOptions) map { credentials =>
      credentials.getAuthenticationType should be("OAuth2")
      secretsMockFile.delete(true)
      dataStoreMockDir.delete(true)
      succeed
    }
  }

  it should "validate" in {
    val secretsMockFile = File.newTemporaryFile("secrets.", ".json").write(GoogleAuthModeSpec.userCredentialsContents)
    val dataStoreMockDir = File.newTemporaryDirectory("dataStore.")
    val userMode = UserMode(
      "user",
      "alice",
      secretsMockFile.pathAsString,
      dataStoreMockDir.pathAsString,
      GoogleConfiguration.GoogleScopes
    )
    val workflowOptions = GoogleAuthModeSpec.emptyOptions
    userMode.validate(workflowOptions)
    secretsMockFile.delete(true)
    dataStoreMockDir.delete(true)
    succeed
  }

  it should "requiresAuthFile" in {
    val secretsMockFile = File.newTemporaryFile("secrets.", ".json").write(GoogleAuthModeSpec.userCredentialsContents)
    val dataStoreMockDir = File.newTemporaryDirectory("dataStore.")
    val userMode = UserMode(
      "user",
      "alice",
      secretsMockFile.pathAsString,
      dataStoreMockDir.pathAsString,
      GoogleConfiguration.GoogleScopes
    )
    userMode.requiresAuthFile should be(false)
    secretsMockFile.delete(true)
    dataStoreMockDir.delete(true)
    succeed
  }

}
