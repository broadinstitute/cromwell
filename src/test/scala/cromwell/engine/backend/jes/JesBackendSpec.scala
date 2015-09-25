package cromwell.engine.backend.jes

import java.net.URL

import cromwell.binding.values.{WdlFile, WdlString}
import cromwell.binding.{Call, CallInputs}
import org.scalatest.mock.MockitoSugar
import org.scalatest.{FlatSpec, Matchers}

class JesBackendSpec extends FlatSpec with Matchers with MockitoSugar {

  "adjustInputPaths" should "map GCS paths and *only* GCS paths to local" in {
    val ignoredCall = mock[Call]
    val stringKey = "abc"
    val stringVal = WdlString("abc")
    val localFileKey = "lf"
    val localFileVal = WdlFile("/blah/abc")
    val gcsFileKey = "gcsf"
    val gcsFileVal = WdlFile("gs://blah/abc")


    val inputs: CallInputs = collection.immutable.HashMap(
      stringKey -> stringVal,
      localFileKey -> localFileVal,
      gcsFileKey -> gcsFileVal
    )

    val mappedInputs: CallInputs  = new JesBackend().adjustInputPaths(ignoredCall, inputs)

    mappedInputs.get(stringKey).get match {
      case WdlString(v) => assert(v.equalsIgnoreCase(stringVal.value))
      case _ => fail("test setup error")
    }

    mappedInputs.get(localFileKey).get match {
      case WdlFile(v) => assert(v.equalsIgnoreCase(localFileVal.value))
      case _ => fail("test setup error")
    }

    mappedInputs.get(gcsFileKey).get match {
      case WdlFile(v) => assert(v.equalsIgnoreCase("/cromwell_root/blah/abc"))
      case _ => fail("test setup error")
    }
  }

  "workflow options existence" should "be verified when in 'RefreshTokenMode'" in {
    val goodOptions = Map("account_name" -> "account", "refresh_token" -> "token")
    val missingToken = Map("account_name" -> "account")
    val missingAccount = Map("refresh_token" -> "token")
    val jesBackend = new JesBackend() {
      override lazy val conf = new JesAttributes(applicationName = "",
      project = "", executionBucket = "", endpointUrl = new URL("http://fakeurl.com"), authMode = RefreshTokenMode, docker = None)
    }

    try {
      jesBackend.assertWorkflowOptions(goodOptions)
    } catch {
      case e: IllegalArgumentException => fail("Correct options validation should not throw an exception.")
      case t: Throwable =>
        t.printStackTrace()
        fail(s"Unexpected exception: ${t.getMessage}")
    }

    the [IllegalArgumentException] thrownBy {
      jesBackend.assertWorkflowOptions(missingToken)
    } should have message s"Missing parameters in workflow options: refresh_token"

    the [IllegalArgumentException] thrownBy {
      jesBackend.assertWorkflowOptions(missingAccount)
    } should have message s"Missing parameters in workflow options: account_name"

    the [IllegalArgumentException] thrownBy {
      jesBackend.assertWorkflowOptions(Map.empty)
    } should have message s"Missing parameters in workflow options: account_name, refresh_token"
  }
}
