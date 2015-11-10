package cromwell.engine.backend.jes

import java.net.URL
import java.nio.file.Paths

import cromwell.binding.CallInputs
import cromwell.binding.types.{WdlArrayType, WdlFileType, WdlMapType, WdlStringType}
import cromwell.binding.values.{WdlArray, WdlFile, WdlMap, WdlString}
import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend.jes.JesBackend.{JesInput, JesOutput}
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

  it should "generate correct JesInputs from a WdlMap" in {
    val inputs = Map(
      "stringToFileMap" -> WdlMap(WdlMapType(WdlStringType, WdlFileType), Map(
        WdlString("stringTofile1") -> WdlFile("gs://path/to/stringTofile1"),
        WdlString("stringTofile2") -> WdlFile("gs://path/to/stringTofile2")
      )),
      "fileToStringMap" -> WdlMap(WdlMapType(WdlFileType, WdlStringType), Map(
        WdlFile("gs://path/to/fileToString1") -> WdlString("fileToString1"),
        WdlFile("gs://path/to/fileToString2") -> WdlString("fileToString2")
      )),
      "fileToFileMap" -> WdlMap(WdlMapType(WdlFileType, WdlFileType), Map(
        WdlFile("gs://path/to/fileToFile1Key") -> WdlFile("gs://path/to/fileToFile1Value"),
        WdlFile("gs://path/to/fileToFile2Key") -> WdlFile("gs://path/to/fileToFile2Value")
      )),
      "stringToString" -> WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
        WdlString("stringToString1") -> WdlString("path/to/stringToString1"),
        WdlString("stringToString2") -> WdlString("path/to/stringToString2")
      ))
    )
    val backendCall = mock[JesBackendCall]
    backendCall.locallyQualifiedInputs returns inputs
    val jesInputs = jesBackend.generateJesInputs(backendCall)
    jesInputs should have size 8
    jesInputs should contain(JesInput("stringToFileMap-0", "gs://path/to/stringTofile1", Paths.get("/cromwell_root/path/to/stringTofile1")))
    jesInputs should contain(JesInput("stringToFileMap-1", "gs://path/to/stringTofile2", Paths.get("/cromwell_root/path/to/stringTofile2")))
    jesInputs should contain(JesInput("fileToStringMap-0", "gs://path/to/fileToString1", Paths.get("/cromwell_root/path/to/fileToString1")))
    jesInputs should contain(JesInput("fileToStringMap-1", "gs://path/to/fileToString2", Paths.get("/cromwell_root/path/to/fileToString2")))
    jesInputs should contain(JesInput("fileToFileMap-0", "gs://path/to/fileToFile1Key", Paths.get("/cromwell_root/path/to/fileToFile1Key")))
    jesInputs should contain(JesInput("fileToFileMap-1", "gs://path/to/fileToFile1Value", Paths.get("/cromwell_root/path/to/fileToFile1Value")))
    jesInputs should contain(JesInput("fileToFileMap-2", "gs://path/to/fileToFile2Key", Paths.get("/cromwell_root/path/to/fileToFile2Key")))
    jesInputs should contain(JesInput("fileToFileMap-3", "gs://path/to/fileToFile2Value", Paths.get("/cromwell_root/path/to/fileToFile2Value")))
  }

  it should "generate correct JesInputs from a WdlArray" in {
    val inputs = Map(
      "fileArray" -> WdlArray(WdlArrayType(WdlFileType), Seq(WdlFile("gs://path/to/file1"), WdlFile("gs://path/to/file2")))
    )
    val backendCall = mock[JesBackendCall]
    backendCall.locallyQualifiedInputs returns inputs
    val jesInputs = jesBackend.generateJesInputs(backendCall)
    jesInputs should have size 2
    jesInputs should contain(JesInput("fileArray-0", "gs://path/to/file1", Paths.get("/cromwell_root/path/to/file1")))
    jesInputs should contain(JesInput("fileArray-1", "gs://path/to/file2", Paths.get("/cromwell_root/path/to/file2")))
  }

  it should "generate correct JesInputs from a WdlFile" in {
    val inputs = Map(
      "file1" -> WdlFile("gs://path/to/file1"),
      "file2" -> WdlFile("gs://path/to/file2")
    )
    val backendCall = mock[JesBackendCall]
    backendCall.locallyQualifiedInputs returns inputs
    val jesInputs = jesBackend.generateJesInputs(backendCall)
    jesInputs should have size 2
    jesInputs should contain(JesInput("file1-0", "gs://path/to/file1", Paths.get("/cromwell_root/path/to/file1")))
    jesInputs should contain(JesInput("file2-0", "gs://path/to/file2", Paths.get("/cromwell_root/path/to/file2")))
  }

  it should "convert local Paths back to corresponding GCS paths in JesOutputs" in {
    val jesOutputs = Seq(
      JesOutput("/cromwell_root/path/to/file1", "gs://path/to/file1", Paths.get("/cromwell_root/path/to/file1")),
      JesOutput("/cromwell_root/path/to/file2", "gs://path/to/file2", Paths.get("/cromwell_root/path/to/file2")),
      JesOutput("/cromwell_root/path/to/file3", "gs://path/to/file3", Paths.get("/cromwell_root/path/to/file3")),
      JesOutput("/cromwell_root/path/to/file4", "gs://path/to/file4", Paths.get("/cromwell_root/path/to/file4")),
      JesOutput("/cromwell_root/path/to/file5", "gs://path/to/file5", Paths.get("/cromwell_root/path/to/file5"))
    )
    val outputValues = Seq(
      WdlFile("/cromwell_root/path/to/file1"),
      WdlArray(WdlArrayType(WdlFileType), Seq(WdlFile("/cromwell_root/path/to/file2"), WdlFile("/cromwell_root/path/to/file3"))),
      WdlMap(WdlMapType(WdlFileType, WdlFileType), Map(
        WdlFile("/cromwell_root/path/to/file4") -> WdlFile("/cromwell_root/path/to/file5")
      ))
    )
    val result = outputValues map jesBackend.wdlValueToGcsPath(jesOutputs)
    result should have size 3
    result should contain(WdlFile("gs://path/to/file1"))
    result should contain(WdlArray(WdlArrayType(WdlFileType),
      Seq(WdlFile("gs://path/to/file2"), WdlFile("gs://path/to/file3")))
    )
    result should contain(WdlMap(WdlMapType(WdlFileType, WdlFileType),
      Map(WdlFile("gs://path/to/file4") -> WdlFile("gs://path/to/file5")))
    )
  }

  it should "create a JesInput for the monitoring script, if specified" in {
    val backendCall = mock[JesBackendCall]
    val wd = mock[WorkflowDescriptor]
    backendCall.workflowDescriptor returns wd

    wd.workflowOptions returns WorkflowOptions.fromJsonString("""{"monitoring_script": "gs://path/to/script"}""").get
    jesBackend.monitoringIO(backendCall) shouldBe Some(JesInput("monitoring", "gs://path/to/script", Paths.get("/cromwell_root/monitoring.sh"),"REFERENCE"))

    wd.workflowOptions returns WorkflowOptions.fromJsonString("""{}""").get
    jesBackend.monitoringIO(backendCall) shouldBe None
  }
}
