package cromwell.engine.backend.jes

import java.net.URL


import cromwell.binding.CallInputs
import cromwell.binding.values.{WdlFile, WdlString}
import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend.jes.authentication._
import cromwell.engine.workflow.{CallKey, WorkflowOptions}
import cromwell.util.EncryptionSpec
import cromwell.util.google.SimpleClientSecrets
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito

import scala.util.Success

class JesBackendSpec extends FlatSpec with Matchers with Mockito {

  val jesBackend = new JesBackend() {
    private val anyString = ""
    private val anyURL: URL = null
    override lazy val jesConf = new JesAttributes(
      applicationName = anyString,
      project = anyString,
      executionBucket = anyString,
      endpointUrl = anyURL,
      authMode = JesAuthMode.fromString("service"),
      dockerCredentials = None,
      googleSecrets = Some(SimpleClientSecrets("myclientId", "myclientSecret"))) {
      override val localizeWithRefreshToken = true
    }
    override lazy val jesConnection = null
  }

  "adjustInputPaths" should "map GCS paths and *only* GCS paths to local" in {
    val ignoredCall = mock[CallKey]
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

    val mappedInputs: CallInputs  = new JesBackend().adjustInputPaths(ignoredCall, inputs, mock[WorkflowDescriptor])

    mappedInputs.get(stringKey).get match {
      case WdlString(v) => assert(v.equalsIgnoreCase(stringVal.value))
      case _ => fail("test setup error")
    }

    mappedInputs.get(localFileKey).get match {
      case wdlFile: WdlFile => assert(wdlFile.value.equalsIgnoreCase(localFileVal.value))
      case _ => fail("test setup error")
    }

    mappedInputs.get(gcsFileKey).get match {
      case wdlFile: WdlFile => assert(wdlFile.value.equalsIgnoreCase("/cromwell_root/blah/abc"))
      case _ => fail("test setup error")
    }
  }

  "workflow options existence" should "be verified when localizing with Refresh Token" in {
    EncryptionSpec.assumeAes256Cbc()

    val goodOptions = WorkflowOptions.fromMap(Map("refresh_token" -> "token")).get
    val missingToken = WorkflowOptions.fromMap(Map.empty).get

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
  }

  it should "create a GcsAuthInformation instance" in {
    val workflowDescriptor = mock[WorkflowDescriptor]
    val mockedWfOptions = mock[WorkflowOptions]
    workflowDescriptor.workflowOptions returns mockedWfOptions
    mockedWfOptions.get("refresh_token") returns Success("myRefreshToken")

    jesBackend.getGcsAuthInformation(workflowDescriptor) shouldBe Some(GcsLocalizing(jesBackend.jesConf.googleSecrets.get, "myRefreshToken"))
  }

}
