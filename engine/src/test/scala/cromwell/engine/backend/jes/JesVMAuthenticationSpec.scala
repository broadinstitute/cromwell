package cromwell.engine.backend.jes

import cromwell.engine.backend.jes.authentication._
import cromwell.engine.io.gcs.SimpleClientSecrets
import org.scalatest.{FlatSpec, Matchers}
import spray.json._


class JesVMAuthenticationSpec extends FlatSpec with Matchers {
  import JesBackend._

  def normalize(str: String) = {
    str.parseJson.prettyPrint
  }

  "generateJson" should "generate the correct json content depending on the available configuration" in {
    val dockerAuth = Option(new JesDockerCredentials("my@docker.account", "mydockertoken"))
    val gcsUserAuth = Option(GcsLocalizing(SimpleClientSecrets("myclientid", "myclientsecret"), "mytoken"))

    generateAuthJson(None, None) shouldBe None

    generateAuthJson(dockerAuth, None) shouldBe
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

    generateAuthJson(None, gcsUserAuth) shouldBe
    Some(normalize("""
           |{
           |    "auths": {
           |        "boto": {
           |            "client_id": "myclientid",
           |            "client_secret": "myclientsecret",
           |            "refresh_token": "mytoken"
           |        }
           |    }
           |}
         """.stripMargin))

    generateAuthJson(dockerAuth, gcsUserAuth) shouldBe
    Some(normalize("""
           |{
           |    "auths": {
           |        "docker": {
           |            "account": "my@docker.account",
           |            "token": "mydockertoken"
           |        },
           |        "boto": {
           |            "client_id": "myclientid",
           |            "client_secret": "myclientsecret",
           |            "refresh_token": "mytoken"
           |        }
           |    }
           |}
         """.stripMargin))

  }

}
