//package cromwell.backend.google.batch.actors
//
//import common.mock.MockSugar
//import cromwell.core.TestKitSuite
//import org.scalatest.flatspec.AnyFlatSpecLike
//import org.scalatest.matchers.should.Matchers
//
//class BatchPollResultMonitorActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with MockSugar {
//
//  behavior of "BatchPollResultMonitorActor"
//
//  it should "set vmCostPerHour if empty" in {
//  }
//
//}

package cromwell.backend.google.batch.actors

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.{TestKit, TestProbe}
import cats.data.Validated.Valid
import common.mock.MockSugar
import cromwell.backend.google.batch.models.GcpBatchRuntimeAttributes
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, RuntimeAttributeDefinition}
import cromwell.core.callcaching.NoDocker
import cromwell.core.{ExecutionEvent, WorkflowOptions}
import cromwell.core.logging.JobLogger
import cromwell.services.cost.{GcpCostLookupRequest, GcpCostLookupResponse, InstantiatedVmInfo}
import cromwell.services.keyvalue.InMemoryKvServiceActor
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import cromwell.backend.google.batch.models.GcpBatchTestConfig._
import wom.graph.CommandCallNode
import cromwell.backend._
import cromwell.backend.google.batch.models._
import cromwell.backend.io.TestWorkflows
import cromwell.backend.standard.pollmonitoring.ProcessThisPollResult
import cromwell.services.metadata.CallMetadataKeys
import cromwell.services.metadata.MetadataService.PutMetadataAction
import org.slf4j.helpers.NOPLogger
import wom.values.WomString

import java.time.{Instant, OffsetDateTime}
import java.time.temporal.ChronoUnit
import scala.concurrent.duration.DurationInt

class BatchPollResultMonitorActorSpec
    extends TestKit(ActorSystem("BatchPollResultMonitorActorSpec"))
    with AnyFlatSpecLike
    with BackendSpec
    with Matchers
    with MockSugar {

  var kvService: ActorRef = system.actorOf(Props(new InMemoryKvServiceActor), "kvService")
  val runtimeAttributesBuilder = GcpBatchRuntimeAttributes.runtimeAttributesBuilder(gcpBatchConfiguration)
  val jobLogger = mock[JobLogger]
  val serviceRegistry = TestProbe()

  val workflowDescriptor = buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld)
  val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
  val jobKey = BackendJobDescriptorKey(call, None, 1)

  val jobDescriptor = BackendJobDescriptor(workflowDescriptor,
                                           jobKey,
                                           runtimeAttributes = Map.empty,
                                           evaluatedTaskInputs = Map.empty,
                                           NoDocker,
                                           None,
                                           Map.empty
  )

  val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"))

  val staticRuntimeAttributeDefinitions: Set[RuntimeAttributeDefinition] =
    GcpBatchRuntimeAttributes.runtimeAttributesBuilder(GcpBatchTestConfig.gcpBatchConfiguration).definitions.toSet

  val defaultedAttributes =
    RuntimeAttributeDefinition.addDefaultsToAttributes(staticRuntimeAttributeDefinitions,
                                                       WorkflowOptions.fromMap(Map.empty).get
    )(
      runtimeAttributes
    )
  val validatedRuntimeAttributes = runtimeAttributesBuilder.build(defaultedAttributes, NOPLogger.NOP_LOGGER)

  val actor = system.actorOf(
    BatchPollResultMonitorActor.props(serviceRegistry.ref,
                                      workflowDescriptor,
                                      jobDescriptor,
                                      validatedRuntimeAttributes,
                                      Some(Gcp),
                                      jobLogger
    )
  )
  val vmInfo = InstantiatedVmInfo("europe-west9", "custom-16-32768", false)


  behavior of "BatchPollResultMonitorActor"

  it should "send a cost lookup request with the correct vm info after receiving a success pollResult" in {

    val terminalPollResult =
      RunStatus.Success(Seq(ExecutionEvent("fakeEvent", OffsetDateTime.now().truncatedTo(ChronoUnit.MILLIS))),
        Some(vmInfo)
      )
    val message = ProcessThisPollResult(terminalPollResult)

    actor ! message

    serviceRegistry.expectMsgPF(1.seconds) { case m: GcpCostLookupRequest =>
      m.vmInfo shouldBe vmInfo
    }
  }

  it should "emit the correct cost metadata after receiving a costLookupResponse" in {

    val costLookupResponse = GcpCostLookupResponse(vmInfo, Valid(BigDecimal(0.1)))

    actor ! costLookupResponse

    serviceRegistry.expectMsgPF(1.seconds) { case m: PutMetadataAction =>
      val event = m.events.head
      m.events.size shouldBe 1
      event.key.key shouldBe CallMetadataKeys.VmCostPerHour
      event.value.get.value shouldBe "0.1"
    }
  }

  it should "emit the correct start time after receiving a running pollResult" in {

    val vmStartTime = OffsetDateTime.now().minus(2, ChronoUnit.HOURS)
    val pollResult = RunStatus.Running(
      Seq(ExecutionEvent(CallMetadataKeys.VmStartTime, vmStartTime)),
      Some(vmInfo)
    )
    val message = ProcessThisPollResult(pollResult)

    actor ! message

    serviceRegistry.expectMsgPF(1.seconds) { case m: PutMetadataAction =>
      val event = m.events.head
      m.events.size shouldBe 1
      event.key.key shouldBe CallMetadataKeys.VmStartTime
      assert(
        Instant
          .parse(event.value.get.value)
          .equals(vmStartTime.toInstant.truncatedTo(ChronoUnit.MILLIS))
      )
    }
  }

  it should "emit the correct end time after receiving a running pollResult" in {

    val vmEndTime = OffsetDateTime.now().minus(2, ChronoUnit.HOURS)
    val pollResult = RunStatus.Running(
      Seq(ExecutionEvent(CallMetadataKeys.VmEndTime, vmEndTime)),
      Some(vmInfo)
    )
    val message = ProcessThisPollResult(pollResult)

    actor ! message

    serviceRegistry.expectMsgPF(1.seconds) { case m: PutMetadataAction =>
      val event = m.events.head
      m.events.size shouldBe 1
      event.key.key shouldBe CallMetadataKeys.VmEndTime
      assert(
        Instant
          .parse(event.value.get.value)
          .equals(vmEndTime.toInstant.truncatedTo(ChronoUnit.MILLIS))
      )
    }
  }
}
