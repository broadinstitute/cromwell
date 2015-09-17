package cromwell.engine.backend.jes

import org.scalatest.mock.MockitoSugar
import org.scalatest.{FlatSpec, Matchers}
import play.api.libs.json._


class GcsAuthSpec extends FlatSpec with Matchers with MockitoSugar {

  def normalize(str: String) = {
    Json.prettyPrint(Json.parse(str))
  }

  "GcsAuthMode" should "parse AuthenticationMode" in {
    GcsAuthMode.fromString("refresh_token") shouldBe RefreshTokenMode
    GcsAuthMode.fromString("service_account") shouldBe ServiceAccountMode
    an [IllegalArgumentException] should be thrownBy GcsAuthMode.fromString("unrecognized_value")
  }

  "generateJson" should "generate the correct json content depending on the available configuration" in {
    val dockerAuth = Option(DockerAuthInformation("my@docker.account", "mydockertoken"))
    val userAuth = Option(GcsUserAuthInformation("my@email.com", "myrefreshtoken"))

    GcsAuth.generateJson(None, None) shouldBe None

    GcsAuth.generateJson(dockerAuth, None) shouldBe
      Some(normalize("""
        |{
        |    "auths": {
        |        "docker": {
        |            "account": "my@docker.account",
        |            "token": "mydockertoken"
        |        }
        |    }
        |}
      """.stripMargin))

    GcsAuth.generateJson(None, userAuth) shouldBe
    Some(normalize("""
           |{
           |    "auths": {
           |        "gcloud": {
           |            "account": "my@email.com",
           |            "token": "myrefreshtoken"
           |        }
           |    }
           |}
         """.stripMargin))

    GcsAuth.generateJson(dockerAuth, userAuth) shouldBe
    Some(normalize("""
           |{
           |    "auths": {
           |        "docker": {
           |            "account": "my@docker.account",
           |            "token": "mydockertoken"
           |        },
           |        "gcloud": {
           |            "account": "my@email.com",
           |            "token": "myrefreshtoken"
           |        }
           |    }
           |}
         """.stripMargin))

  }

}
