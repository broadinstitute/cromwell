package cromwell.core

import cromwell.util.EncryptionSpec
import org.scalatest.{Matchers, WordSpecLike}
import spray.json._

import scala.util.{Failure, Success}

class WorkflowOptionsSpec extends Matchers with WordSpecLike {
  val workflowOptionsJson =
    """{
      |  "key": "value",
      |  "refresh_token": "it's a secret"
      |}
    """.stripMargin.parseJson.asInstanceOf[JsObject]

  "WorkflowOptions" should {
    "parse workflow options properly" in {
      EncryptionSpec.assumeAes256Cbc()

      WorkflowOptions.fromJsonObject(workflowOptionsJson) match {
        case Success(options) =>
          options.get("key") shouldEqual Success("value")
          options.get("refresh_token") shouldEqual Success("it's a secret")
          options.get("bad_key") shouldBe a [Failure[_]]

          val encryptedOptions = options.asPrettyJson.parseJson.asInstanceOf[JsObject]
          val refreshTokenEncrypted = encryptedOptions.fields.getOrElse("refresh_token", fail("refresh_token key expected"))
          refreshTokenEncrypted shouldBe a [JsObject]
          refreshTokenEncrypted.asInstanceOf[JsObject].fields.keys shouldEqual Set("iv", "ciphertext")

          options.clearEncryptedValues shouldEqual """{
                                                |  "key": "value",
                                                |  "refresh_token": "cleared"
                                                |}""".stripMargin
        case _ => fail("Expecting workflow options to be parseable")
      }
    }
  }
}