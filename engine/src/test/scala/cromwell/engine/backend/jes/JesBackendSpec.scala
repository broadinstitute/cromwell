package cromwell.engine.backend.jes

import cromwell.CromwellTestkitSpec
import cromwell.backend.impl.jes.io.{DiskType, JesWorkingDisk}
import cromwell.core.WorkflowOptions
import cromwell.engine.backend.EnhancedWorkflowOptions._
import cromwell.engine.workflow.WorkflowDescriptorBuilder
import cromwell.filesystems.gcs._
import cromwell.util.EncryptionSpec
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito

class JesBackendSpec extends FlatSpec with Matchers with Mockito with BeforeAndAfterAll with WorkflowDescriptorBuilder {
  val testWorkflowManagerSystem = new CromwellTestkitSpec.TestWorkflowManagerSystem()
  override implicit val actorSystem = testWorkflowManagerSystem.actorSystem
  val workingDisk = JesWorkingDisk(DiskType.SSD, 200)

  override protected def afterAll() = {
    testWorkflowManagerSystem.shutdownTestActorSystem()
    super.afterAll()
  }

  val refreshToken = RefreshTokenMode(name = "bar", clientId = "secret-id", clientSecret = "secret-secret")

  "workflow options existence" should "be verified when localizing with Refresh Token" in {
    EncryptionSpec.assumeAes256Cbc()

    val goodOptions = WorkflowOptions.fromMap(Map("refresh_token" -> "token")).get
    refreshToken.assertWorkflowOptions(goodOptions.toGoogleAuthOptions)

    val badOptions = WorkflowOptions.fromMap(Map("fresh_tokin" -> "broken")).get
    val noOptions = WorkflowOptions.fromMap(Map.empty[String, String]).get

    List(badOptions, noOptions).foreach { option =>
      the [IllegalArgumentException] thrownBy {
        refreshToken.assertWorkflowOptions(option.toGoogleAuthOptions)
      } should have message s"Missing parameters in workflow options: refresh_token"
    }
  }
}
