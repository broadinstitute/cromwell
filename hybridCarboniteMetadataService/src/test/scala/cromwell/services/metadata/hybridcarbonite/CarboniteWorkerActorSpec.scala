package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorRef
import akka.testkit.{TestActorRef, TestProbe}
import com.typesafe.config.ConfigFactory
import common.validation.Validation._
import cromwell.core.{TestKitSuite, WorkflowAborted, WorkflowFailed, WorkflowId, WorkflowSucceeded}
import cromwell.services.metadata.MetadataArchiveStatus.Unarchived
import cromwell.services.metadata.MetadataService.{QueryForWorkflowsMatchingParameters, QueryMetadata, WorkflowQueryResponse, WorkflowQueryResult, WorkflowQuerySuccess}
import cromwell.services.metadata.WorkflowQueryKey._
import cromwell.services.metadata.hybridcarbonite.CarboniteWorkerActor.CarboniteWorkflowComplete
import cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActor.FreezeMetadata
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._

class CarboniteWorkerActorSpec extends TestKitSuite("CarboniteWorkerActorSpec") with FlatSpecLike with Matchers {

  val serviceRegistryActor = TestProbe()
  val ioActor = TestProbe()
  val carboniteFreezerActor = TestProbe()

  val carboniterConfig = HybridCarboniteConfig.parseConfig(ConfigFactory.parseString(
    """{
      |   enabled = true
      |   bucket = "carbonite-test-bucket"
      |   filesystems {
      |     gcs {
      |       # A reference to the auth to use for storing and retrieving metadata:
      |       auth = "application-default"
      |     }
      |   }
      |}""".stripMargin
  )).unsafe("Make config file")

  val findWorkflowToCarboniteQuery: Seq[(String, String)] = Seq(
    IncludeSubworkflows.name -> "false",
    Status.name -> WorkflowSucceeded.toString,
    Status.name -> WorkflowFailed.toString,
    Status.name -> WorkflowAborted.toString,
    MetadataArchiveStatus.name -> Unarchived.toString,
    Page.name -> "1",
    PageSize.name -> "1"
  )

  it should "carbonite workflow at intervals" in {
    val carboniteWorkerActor = TestActorRef(new MockCarboniteWorkerActor(carboniterConfig, serviceRegistryActor.ref, ioActor.ref, carboniteFreezerActor.ref))

    val workflowToCarbonite = "04c93860-ea0a-11e9-81b4-2a2ae2dbcce4"
    val queryMeta = QueryMetadata(Option(1), Option(1), Option(1))
    val queryResult = WorkflowQueryResult(workflowToCarbonite, None, None, None, None, None, None, None, None, Unarchived)
    val queryResponse = WorkflowQueryResponse(Seq(queryResult), 1)
    val querySuccessResponse = WorkflowQuerySuccess(queryResponse, Option(queryMeta))

    serviceRegistryActor.expectMsgPF(10.seconds) {
      case QueryForWorkflowsMatchingParameters(queryParams) => queryParams shouldBe findWorkflowToCarboniteQuery
    }

    serviceRegistryActor.send(carboniteWorkerActor, querySuccessResponse)

    carboniteFreezerActor.expectMsgPF(10.seconds) {
      case FreezeMetadata(id) =>
        id.toString shouldBe workflowToCarbonite
    }

    carboniteFreezerActor.send(carboniteWorkerActor, CarboniteWorkflowComplete(WorkflowId.fromString(workflowToCarbonite)))

    /*
      after carbonation of a workflow is complete, the CarboniteWorkerActor should schedule to ask ServiceRegistry for
      the next workflow to carbonite
     */
    serviceRegistryActor.expectMsgPF(10.seconds) {
      case QueryForWorkflowsMatchingParameters(queryParams) => queryParams shouldBe findWorkflowToCarboniteQuery
    }
  }
}



class MockCarboniteWorkerActor(carboniterConfig: HybridCarboniteConfig,
                               serviceRegistryActor: ActorRef,
                               ioActor: ActorRef,
                               override val carboniteFreezerActor: ActorRef)
  extends CarboniteWorkerActor(carboniterConfig, serviceRegistryActor, ioActor) { }
