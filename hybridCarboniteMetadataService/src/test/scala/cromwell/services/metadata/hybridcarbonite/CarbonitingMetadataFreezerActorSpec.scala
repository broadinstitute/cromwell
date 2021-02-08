package cromwell.services.metadata.hybridcarbonite

import akka.actor.{ActorRef, FSM}
import akka.testkit.{TestFSMRef, TestProbe}
import com.typesafe.config.ConfigFactory
import common.assertion.ManyTimes._
import common.validation.Validation._
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise
import cromwell.core.io.IoWriteCommand
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataService.GetMetadataAction
import cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActor._
import cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActorSpec.TestableCarbonitingMetadataFreezerActor
import cromwell.services.{FailedMetadataJsonResponse, SuccessfulMetadataJsonResponse}
import org.scalatest.concurrent.Eventually._
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import spray.json._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}

class CarbonitingMetadataFreezerActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers {

  implicit val ec: ExecutionContext = system.dispatcher

  private val carboniterConfigString =
    """
      |bucket = "carbonite-test-bucket"
      |filesystems {
      |  gcs {
      |    # A reference to the auth to use for storing and retrieving metadata:
      |    auth = "application-default"
      |  }
      |}
      |metadata-freezing {
      |  initial-interval: 5 seconds
      |}
      |[BUCKETREADLIMIT]
      |""".stripMargin

  val carboniterConfigDefaultLimit: HybridCarboniteConfig =
    HybridCarboniteConfig
      .parseConfig(
        ConfigFactory.parseString(carboniterConfigString.replace("[BUCKETREADLIMIT]", ""))
      ).unsafe("Make config file")

  // any metadata will be too big with this config
  val carboniterConfigZeroBytesLimit: HybridCarboniteConfig =
    HybridCarboniteConfig
      .parseConfig(
        ConfigFactory.parseString(carboniterConfigString.replace("[BUCKETREADLIMIT]", "bucket-read-limit-bytes = 0"))
      ).unsafe("Make config file")

  val serviceRegistryActor: TestProbe = TestProbe("serviceRegistryActor")
  val ioActor: TestProbe = TestProbe("ioActor")
  val carboniteWorkerActor: TestProbe = TestProbe("carboniteWorkerActor")

  it should "follow the expected golden-path lifecycle" in {
    val actor = TestFSMRef(new TestableCarbonitingMetadataFreezerActor(carboniterConfigDefaultLimit, carboniteWorkerActor.ref, serviceRegistryActor.ref, ioActor.ref))
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

    serviceRegistryActor.send(actor, SuccessfulMetadataJsonResponse(null, jsonValue.asJsObject))

    var ioCommandPromise: Promise[Any] = null
    ioActor.expectMsgPF(10.seconds) {
      case command @ IoCommandWithPromise(ioCommand: IoWriteCommand, _) =>
        ioCommand.file.pathAsString should be(HybridCarboniteConfig.pathForWorkflow(workflowIdToFreeze, "carbonite-test-bucket"))
        ioCommand.content should be(jsonValue.prettyPrint)
        ioCommand.compressPayload should be(true)
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
    eventually {
      actor.stateName should be(Pending)
    }
  }

  it should "mark workflow as TooLargeToArchive in case if resulting json size exceeds the configured bucket-read-limit-bytes value" in {
    val actor = TestFSMRef(new TestableCarbonitingMetadataFreezerActor(carboniterConfigZeroBytesLimit, carboniteWorkerActor.ref, serviceRegistryActor.ref, ioActor.ref))
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

    serviceRegistryActor.send(actor, SuccessfulMetadataJsonResponse(null, jsonValue.asJsObject))

    eventually {
      actor.underlyingActor.updateArchiveStatusCall should be((workflowIdToFreeze, MetadataArchiveStatus.TooLargeToArchive))
      actor.stateName should be(UpdatingDatabase)
    }

    // Finally, when the database responds appropriately, the actor should shut itself down:
    watch(actor)
    actor.underlyingActor.updateArchiveStatusPromise.success(1)
    eventually {
      actor.stateName should be(Pending)
    }
  }

  it should "correctly handle failure received from IOActor" in {
    val actor = TestFSMRef(new TestableCarbonitingMetadataFreezerActor(carboniterConfigDefaultLimit, carboniteWorkerActor.ref, serviceRegistryActor.ref, ioActor.ref))
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

    serviceRegistryActor.send(actor, SuccessfulMetadataJsonResponse(null, jsonValue.asJsObject))

    var ioCommandPromise: Promise[Any] = null
    ioActor.expectMsgPF(10.seconds) {
      case command @ IoCommandWithPromise(ioCommand: IoWriteCommand, _) =>
        ioCommand.file.pathAsString should be(HybridCarboniteConfig.pathForWorkflow(workflowIdToFreeze, "carbonite-test-bucket"))
        ioCommand.content should be(jsonValue.prettyPrint)
        ioCommand.compressPayload should be(true)
        ioCommandPromise = command.promise
    }
    eventually {
      actor.stateName should be(Freezing)
    }

    // Simulates the IoActor failing to complete the IO command successfully:
    ioCommandPromise.failure(new RuntimeException("TEST EXCEPTION: Cannot write metadata to GCS bucket"))
    eventually {
      actor.underlyingActor.updateArchiveStatusCall should be((workflowIdToFreeze, MetadataArchiveStatus.ArchiveFailed))
      actor.stateName should be(UpdatingDatabase)
    }

    // Finally, when the database responds appropriately, the actor should shut itself down:
    watch(actor)
    actor.underlyingActor.updateArchiveStatusPromise.success(1)
    eventually {
      actor.stateName should be(Pending)
    }
  }

  it should "correctly handle failure received on attempt to fetch metadata json for Carboniting" in {
    val actor = TestFSMRef(new TestableCarbonitingMetadataFreezerActor(carboniterConfigDefaultLimit, carboniteWorkerActor.ref, serviceRegistryActor.ref, ioActor.ref))
    val workflowIdToFreeze = WorkflowId.randomId()

    actor ! FreezeMetadata(workflowIdToFreeze)

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case gma: GetMetadataAction => gma.workflowId should be(workflowIdToFreeze)
    }
    eventually {
      actor.stateName should be(Fetching)
    }

    // Freezer actor will receive `FailedMetadataJsonResponse` from service registry actor on any failure, including DB read timeout
    serviceRegistryActor.send(actor, FailedMetadataJsonResponse(null, new RuntimeException("TEST EXCEPTION: Cannot read classic metadata from DB")))

    eventually {
      actor.underlyingActor.updateArchiveStatusCall should be((workflowIdToFreeze, MetadataArchiveStatus.ArchiveFailed))
      actor.stateName should be(UpdatingDatabase)
    }

    // Finally, when the database responds appropriately, the actor should shut itself down:
    watch(actor)
    actor.underlyingActor.updateArchiveStatusPromise.success(1)
    eventually {
      actor.stateName should be(Pending)
    }
  }

  it should "when status update in DB fails retry it until success" in {
    val actor = TestFSMRef(new TestableCarbonitingMetadataFreezerActor(carboniterConfigDefaultLimit, carboniteWorkerActor.ref, serviceRegistryActor.ref, ioActor.ref))
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

    serviceRegistryActor.send(actor, SuccessfulMetadataJsonResponse(null, jsonValue.asJsObject))

    var ioCommandPromise: Promise[Any] = null
    ioActor.expectMsgPF(10.seconds) {
      case command @ IoCommandWithPromise(ioCommand: IoWriteCommand, _) =>
        ioCommand.file.pathAsString should be(HybridCarboniteConfig.pathForWorkflow(workflowIdToFreeze, "carbonite-test-bucket"))
        ioCommand.content should be(jsonValue.prettyPrint)
        ioCommand.compressPayload should be(true)
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

    // When the database update fails, the actor should retry update endlessly until success:
    watch(actor)
    10.times {
      actor.underlyingActor.updateArchiveStatusPromise.failure(new RuntimeException("TEST EXCEPTION: Cannot update status in DB"))
      eventually {
        actor.stateName should be(UpdatingDatabase)
      }
      ()
    }
    actor.underlyingActor.updateArchiveStatusPromise.success(1)
    eventually {
      actor.stateName should be(Pending)
    }
  }

}

object CarbonitingMetadataFreezerActorSpec {

  class TestableCarbonitingMetadataFreezerActor(config: HybridCarboniteConfig, carboniteWorkerActor: ActorRef, serviceRegistry: ActorRef, ioActor: ActorRef)
    extends CarbonitingMetadataFreezerActor(config.freezingConfig.asInstanceOf[ActiveMetadataFreezingConfig], config, carboniteWorkerActor, serviceRegistry, ioActor) {

    var updateArchiveStatusCall: (WorkflowId, MetadataArchiveStatus) = (null, null)
    var updateArchiveStatusPromise: Promise[Int] = Promise[Int]()

    override def updateMetadataArchiveStatus(workflowId: WorkflowId, newStatus: MetadataArchiveStatus): Future[Int] = {
      updateArchiveStatusCall = (workflowId, newStatus)
      updateArchiveStatusPromise.future
    }

    override def scheduleDatabaseUpdateAndAwaitResult(workflowId: WorkflowId,
                                                      newStatus: MetadataArchiveStatus,
                                                      delay: Option[FiniteDuration]): FSM.State[CarbonitingMetadataFreezingState, CarbonitingMetadataFreezingData] = {
      updateArchiveStatusPromise = Promise[Int]()
      super.scheduleDatabaseUpdateAndAwaitResult(workflowId, newStatus, None)
    }
  }

}
