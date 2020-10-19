package cromwell.cloudsupport.gcp.auth

import java.io.FileNotFoundException

import better.files.File
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class UserModeSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "UserMode"

  it should "fail to generate a bad credential" in {
    val secretsMockFile = File.newTemporaryFile("secrets.", ".json").write(GoogleAuthModeSpec.userCredentialsContents)
    val userMode = UserMode("user", secretsMockFile.pathAsString)
    val exception = intercept[RuntimeException](userMode.credentials())
    exception.getMessage should startWith("Google credentials are invalid: ")
    secretsMockFile.delete(swallowIOExceptions = true)
  }

  it should "fail to generate a bad credential from a secrets json" in {
    val secretsMockFile = File.newTemporaryFile("secrets.", ".json").delete()
    val userMode = UserMode("user", secretsMockFile.pathAsString)
    val exception = intercept[FileNotFoundException](userMode.credentials())
    exception.getMessage should fullyMatch regex "File .*/secrets..*.json does not exist or is not readable"
  }

  it should "generate a non-validated credential" in {
    val secretsMockFile = File.newTemporaryFile("secrets.", ".json").write(GoogleAuthModeSpec.userCredentialsContents)
    val userMode = UserMode("user", secretsMockFile.pathAsString)
    userMode.credentialsValidation = GoogleAuthMode.NoCredentialsValidation
    val credentials = userMode.credentials()
    credentials.getAuthenticationType should be("OAuth2")
    secretsMockFile.delete(swallowIOExceptions = true)
  }

  it should "requiresAuthFile" in {
    val secretsMockFile = File.newTemporaryFile("secrets.", ".json").write(GoogleAuthModeSpec.userCredentialsContents)
    val userMode = UserMode("user", secretsMockFile.pathAsString)
    userMode.requiresAuthFile should be(false)
    secretsMockFile.delete(swallowIOExceptions = true)
  }

}
