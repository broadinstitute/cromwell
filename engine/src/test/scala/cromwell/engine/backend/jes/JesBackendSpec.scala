package cromwell.engine.backend.jes

import java.net.URL
import java.nio.file.{FileSystems, Paths}
import java.util.UUID

import akka.actor.ActorSystem
import com.google.api.client.testing.http.{HttpTesting, MockHttpTransport, MockLowLevelHttpRequest, MockLowLevelHttpResponse}
import cromwell.CromwellTestkitSpec
import cromwell.backend.PreemptedException
import cromwell.backend.impl.jes.io.{DiskType, JesWorkingDisk}
import cromwell.core.{OldCallContext, OldWorkflowContext, WorkflowId, WorkflowOptions}
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.filesystems.gcs._
import cromwell.util.{EncryptionSpec, SampleWdl}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.slf4j.{Logger, LoggerFactory}
import org.specs2.mock.Mockito
import wdl4s.types.{WdlArrayType, WdlFileType, WdlMapType, WdlStringType}
import wdl4s.values._
import wdl4s.{Call, CallInputs, Task}

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.util.{Success, Try}

class JesBackendSpec extends FlatSpec with Matchers with Mockito with BeforeAndAfterAll with WorkflowDescriptorBuilder {
  val testWorkflowManagerSystem = new CromwellTestkitSpec.TestWorkflowManagerSystem()
  override implicit val actorSystem = testWorkflowManagerSystem.actorSystem
  val workingDisk = JesWorkingDisk(DiskType.SSD, 200)

  override protected def afterAll() = {
    testWorkflowManagerSystem.shutdownTestActorSystem()
    super.afterAll()
  }

  val clientSecrets = RefreshTokenMode(name = "bar", clientId = "secret-id", clientSecret = "secret-secret")
//  val jesBackend = new OldStyleJesBackend(CromwellTestkitSpec.JesBackendConfigEntry, actorSystem) {
//    private val anyString = ""
//    private val anyURL: URL = null
//    override lazy val jesAttributes = new JesAttributes(
//      project = anyString,
//      executionBucket = anyString,
//      endpointUrl = anyURL,
//      maxPollingInterval = 600,
//      genomicsAuth = ApplicationDefaultMode(name = "foo"),
//      gcsFilesystemAuth = clientSecrets)
//  }

  //"workflow options existence" should "be verified when localizing with Refresh Token" in {
  ignore should "be verified when localizing with Refresh Token" in {
//    EncryptionSpec.assumeAes256Cbc()
//
//    val goodOptions = WorkflowOptions.fromMap(Map("refresh_token" -> "token")).get
//
//    try {
//      jesBackend.assertWorkflowOptions(goodOptions)
//    } catch {
//      case e: IllegalArgumentException => fail("Correct options validation should not throw an exception.")
//      case t: Throwable =>
//        t.printStackTrace()
//        fail(s"Unexpected exception: ${t.getMessage}")
//    }
//
//    val missingToken = WorkflowOptions.fromMap(Map.empty).get
//    the [IllegalArgumentException] thrownBy {
//      jesBackend.assertWorkflowOptions(missingToken)
//    } should have message s"Missing parameters in workflow options: refresh_token"
  }

  ignore should "create a GcsAuthInformation instance" in {
//    val workflowDescriptor = mock[OldStyleWorkflowDescriptor]
//    val mockedWfOptions = mock[WorkflowOptions]
//    workflowDescriptor.workflowOptions returns mockedWfOptions
//    mockedWfOptions.get("refresh_token") returns Success("myRefreshToken")
//
//    jesBackend.refreshTokenAuth(workflowDescriptor) shouldBe Some(GcsLocalizing(clientSecrets, "myRefreshToken"))
  }
}
