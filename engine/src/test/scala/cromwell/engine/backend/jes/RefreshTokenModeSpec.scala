package cromwell.engine.backend.jes

import cromwell.backend.BackendSpec._
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.engine.backend.EnhancedWorkflowOptions._
import cromwell.filesystems.gcs._
import cromwell.util.EncryptionSpec
import org.scalatest.FlatSpecLike

class RefreshTokenModeSpec extends TestKitSuite("RefreshTokenModeSpec") with FlatSpecLike {

  val refreshToken = RefreshTokenMode(name = "bar", clientId = "secret-id", clientSecret = "secret-secret")
  val mockToken = "token"

  "workflow options existence" should "be verified when localizing with Refresh Token" in {
    EncryptionSpec.assumeAes256Cbc()

    val badOptions = WorkflowOptions.fromMap(Map("fresh_tokin" -> "broken")).get
    val noOptions = WorkflowOptions.fromMap(Map.empty[String, String]).get

    List(badOptions, noOptions).foreach { option =>
      the [IllegalArgumentException] thrownBy {
        refreshToken.assertWorkflowOptions(option.toGoogleAuthOptions)
      } should have message s"Missing parameters in workflow options: refresh_token"
    }
  }

  "refresh token value" should "be decrypted when asserting workflow options" in {
    EncryptionSpec.assumeAes256Cbc()

    val goodOptions = WorkflowOptions.fromMap(Map("refresh_token" -> mockToken)).get
    refreshToken.assertWorkflowOptions(goodOptions.toGoogleAuthOptions)
    goodOptions.toGoogleAuthOptions.get("refresh_token").get shouldBe mockToken
  }
}
