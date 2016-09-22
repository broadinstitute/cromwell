package cromwell.filesystems.gcs

import org.scalatest.{FlatSpec, Matchers}

class RefreshTokenModeSpec extends FlatSpec with Matchers {

  val refreshToken = RefreshTokenMode(name = "bar", clientId = "secret-id", clientSecret = "secret-secret")

  behavior of "RefreshTokenMode"

  it should "assert good workflow options" in {
    val goodOptions = GoogleOptionsMap(Map("refresh_token" -> "token"))
    refreshToken.assertWorkflowOptions(goodOptions)
  }

  it should "fail to assert bad workflow options" in {
    val badOptions = GoogleOptionsMap(Map("fresh_tokin" -> "broken"))
    val noOptions = GoogleOptionsMap(Map.empty[String, String])

    List(badOptions, noOptions).foreach { option =>
      the[IllegalArgumentException] thrownBy {
        refreshToken.assertWorkflowOptions(option)
      } should have message s"Missing parameters in workflow options: refresh_token"
    }
  }
}
