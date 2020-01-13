package cromwell.services.metadata.hybridcarbonite

import akka.actor.ActorRef
import akka.testkit.{EventFilter, TestActorRef, TestProbe}
import com.typesafe.config.ConfigFactory
import common.validation.Validation._
import common.assertion.ManyTimes._
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.services.metadata.MetadataArchiveStatus.{Archived, Unarchived}
import cromwell.services.metadata.MetadataService.{DeleteMetadataAction, DeleteMetadataFailedResponse, DeleteMetadataSuccessfulResponse, QueryForWorkflowsMatchingParameters, QueryMetadata, WorkflowQueryResponse, WorkflowQueryResult, WorkflowQuerySuccess}
import cromwell.services.metadata.hybridcarbonite.CarboniteWorkerActor.CarboniteWorkflowComplete
import cromwell.services.metadata.hybridcarbonite.CarbonitingMetadataFreezerActor.FreezeMetadata
import org.scalatest.{FlatSpecLike, Matchers}

import scala.Option
import scala.concurrent.duration._

class CarboniteWorkerActorSpec extends TestKitSuite("CarboniteWorkerActorSpec") with FlatSpecLike with Matchers {

  val serviceRegistryActor = TestProbe()
  val ioActor = TestProbe()
  val carboniteFreezerActor = TestProbe()
  val deleteMetadataActor = TestProbe()

  val fasterBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 10.millis,
    maxInterval = 100.millis,
    multiplier = 1,
    randomizationFactor = 0.0
  )

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

  val carboniteWorkerActor = TestActorRef(new MockCarboniteWorkerActor(
    carboniterConfig,
    serviceRegistryActor.ref,
    ioActor.ref,
    carboniteFreezerActor.ref,
    fasterBackOff
  ))

  val workflowToCarbonite = "04c93860-ea0a-11e9-81b4-2a2ae2dbcce4"
  val queryMeta = QueryMetadata(Option(1), Option(1), Option(1))
  val queryResult = WorkflowQueryResult(workflowToCarbonite, None, None, None, None, None, None, None, None, Unarchived)
  val queryResponse = WorkflowQueryResponse(Seq(queryResult), 1)
  val querySuccessResponse = WorkflowQuerySuccess(queryResponse, Option(queryMeta))

  it should "carbonite workflow at intervals and delete data from DB after carboniting finished" in {
    10.times {
      // We might get noise from instrumentation. We can ignore that, but we expect the query to come through eventually:
      val expectedQueryParams = CarboniteWorkerActor.buildQueryParametersForWorkflowToCarboniteQuery(carboniterConfig.minimumSummaryEntryId)
      serviceRegistryActor.fishForSpecificMessage(10.seconds) {
        case QueryForWorkflowsMatchingParameters(`expectedQueryParams`) => true
      }

      serviceRegistryActor.send(carboniteWorkerActor, querySuccessResponse)

      carboniteFreezerActor.expectMsg(FreezeMetadata(WorkflowId.fromString(workflowToCarbonite)))

      carboniteFreezerActor.send(carboniteWorkerActor, CarboniteWorkflowComplete(WorkflowId.fromString(workflowToCarbonite), Archived))

      serviceRegistryActor.fishForSpecificMessage(10.second) {
        case DeleteMetadataAction(_, _, _) => true
      }

      EventFilter.info(message = s"Completed deleting metadata from database for carbonited workflow: $workflowToCarbonite", occurrences = 1) intercept {
        deleteMetadataActor.send(carboniteWorkerActor, DeleteMetadataSuccessfulResponse(WorkflowId.fromString(workflowToCarbonite)))
      }
    }
  }

  it should "carbonite workflow at intervals even if metadata deletion after carboniting fails" in {
    10.times {
      // We might get noise from instrumentation. We can ignore that, but we expect the query to come through eventually:
      serviceRegistryActor.fishForSpecificMessage(10.seconds) {
        case QueryForWorkflowsMatchingParameters(params) => params == CarboniteWorkerActor.buildQueryParametersForWorkflowToCarboniteQuery(Option(0))
      }

      serviceRegistryActor.send(carboniteWorkerActor, querySuccessResponse)

      carboniteFreezerActor.expectMsg(FreezeMetadata(WorkflowId.fromString(workflowToCarbonite)))

      carboniteFreezerActor.send(carboniteWorkerActor, CarboniteWorkflowComplete(WorkflowId.fromString(workflowToCarbonite), Archived))

      serviceRegistryActor.fishForSpecificMessage(10.second) {
        case DeleteMetadataAction(_, _, _) => true
      }

      EventFilter.error(message = s"All attempts to delete metadata from database for carbonited workflow $workflowToCarbonite failed.", occurrences = 1) intercept {
        deleteMetadataActor.send(carboniteWorkerActor,
          DeleteMetadataFailedResponse(
            WorkflowId.fromString(workflowToCarbonite),
            new RuntimeException("Test exception")
          )
        )
      }
    }
  }
}

class MockCarboniteWorkerActor(carboniterConfig: HybridCarboniteConfig,
                               serviceRegistryActor: ActorRef,
                               ioActor: ActorRef,
                               override val carboniteFreezerActor: ActorRef,
                               override val backOff: SimpleExponentialBackoff)
  extends CarboniteWorkerActor(carboniterConfig, serviceRegistryActor, ioActor) { }
