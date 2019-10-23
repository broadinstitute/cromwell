package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorRef
import akka.testkit.{TestFSMRef, TestProbe}
import com.typesafe.config.ConfigFactory
import common.validation.Validation._
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise
import cromwell.core.io.IoWriteCommand
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataService.GetMetadataAction
import cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActor.{Fetching, FreezeMetadata, Freezing, Pending, UpdatingDatabase}
import cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActorSpec.TestableCarbonitingMetadataFreezerActor
import cromwell.services.BuiltMetadataResponse
import org.scalatest.{FlatSpecLike, Matchers}
import spray.json._
import org.scalatest.concurrent.Eventually._

import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.concurrent.duration._

class CarbonitingMetadataFreezerActorSpec extends TestKitSuite("CarbonitedMetadataThawingActorSpec") with FlatSpecLike with Matchers {

  implicit val ec: ExecutionContext = system.dispatcher

  val carboniterConfig = HybridCarboniteConfig.parseConfig(ConfigFactory.parseString(
    """enabled = true
      |bucket = "carbonite-test-bucket"
      |filesystems {
      |  gcs {
      |    # A reference to the auth to use for storing and retrieving metadata:
      |    auth = "application-default"
      |  }
      |}""".stripMargin)).unsafe("Make config file")

  val serviceRegistryActor = TestProbe()
  val ioActor = TestProbe()
  val carboniteWorkerActor = TestProbe()

  it should "follow the expected golden-path lifecycle" in {
    val actor = TestFSMRef(new TestableCarbonitingMetadataFreezerActor(carboniterConfig, carboniteWorkerActor.ref, serviceRegistryActor.ref, ioActor.ref))
    val workflowIdToFreeze = WorkflowId.randomId()

    actor ! FreezeMetadata(workflowIdToFreeze)

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case gma: GetMetadataAction => gma.workflowId should be(workflowIdToFreeze)
    }
    eventually {
      actor.stateName should be(Fetching)
    }

    val jsonValue =
      s"""{
         |  "id": "$workflowIdToFreeze",
         |  "status": "Successful"
         |}""".stripMargin.parseJson

    serviceRegistryActor.send(actor, BuiltMetadataResponse(null, jsonValue.asJsObject))

    var ioCommandPromise: Promise[Any] = null
    ioActor.expectMsgPF(10.seconds) {
      case command @ IoCommandWithPromise(ioCommand: IoWriteCommand, _) =>
        ioCommand.file.pathAsString should be(HybridCarboniteConfig.pathForWorkflow(workflowIdToFreeze, "carbonite-test-bucket"))
        ioCommand.content should be(jsonValue.prettyPrint)
        ioCommandPromise = command.promise
    }
    eventually {
      actor.stateName should be(Freezing)
    }

    // Simulates the IoActor completing the IO command successfully:
    ioCommandPromise.success(())
    eventually {
      actor.underlyingActor.updateArchiveStatusCall should be((workflowIdToFreeze, MetadataArchiveStatus.Archived))
      actor.stateName should be(UpdatingDatabase)
    }

    // Finally, when the database responds appropriately, the actor should shut itself down:
    watch(actor)
    actor.underlyingActor.updateArchiveStatusPromise.success(1)
    actor.stateName should be(Pending)
  }

}

object CarbonitingMetadataFreezerActorSpec {

  class TestableCarbonitingMetadataFreezerActor(config: HybridCarboniteConfig, carboniteWorkerActor: ActorRef, serviceRegistry: ActorRef, ioActor: ActorRef)
    extends CarbonitingMetadataFreezerActor(config, carboniteWorkerActor, serviceRegistry, ioActor) {

    var updateArchiveStatusCall: (WorkflowId, MetadataArchiveStatus) = (null, null)
    val updateArchiveStatusPromise = Promise[Int]()

    override def updateMetadataArchiveStatus(workflowId: WorkflowId, newStatus: MetadataArchiveStatus): Future[Int] = {
      updateArchiveStatusCall = (workflowId, newStatus)
      updateArchiveStatusPromise.future
    }
  }

}
