package cromwell.core

import cromwell.util.EncryptionSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import spray.json._

import scala.util.{Failure, Success}

class WorkflowOptionsSpec extends Matchers with AnyWordSpecLike {
  val workflowOptionsJson =
    """{
      |  "key": "value"
      |}
    """.stripMargin.parseJson.asInstanceOf[JsObject]

  "WorkflowOptions" should {
    "parse workflow options properly" in {
      EncryptionSpec.assumeAes256Cbc()

      WorkflowOptions.fromJsonObject(workflowOptionsJson) match {
        case Success(options) =>
          options.get("key") shouldEqual Success("value")
          options.get("bad_key") shouldBe a[Failure[_]]

          options.clearEncryptedValues.asPrettyJson shouldEqual """{
                                                                  |  "key": "value"
                                                                  |}""".stripMargin
        case _ => fail("Expecting workflow options to be parseable")
      }
    }
  }
}
