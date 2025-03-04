package cromwell.backend.google.batch.api.request

import akka.actor.ActorSystem
import akka.testkit.TestKit
import com.google.cloud.batch.v1.AllocationPolicy.{
  InstancePolicy,
  InstancePolicyOrTemplate,
  LocationPolicy,
  ProvisioningModel
}
import com.google.cloud.batch.v1.JobStatus.State
import com.google.cloud.batch.v1._
import com.google.protobuf.Timestamp
import common.mock.MockSugar
import cromwell.backend.google.batch.api.BatchApiResponse
import cromwell.backend.google.batch.models.RunStatus
import org.mockito.ArgumentMatchers.any
import org.mockito.Mockito.doReturn
import org.scalatest.PrivateMethodTester
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.jdk.CollectionConverters._

class BatchRequestExecutorSpec
    extends TestKit(ActorSystem("BatchRequestExecutorSpec"))
    with AnyFlatSpecLike
    with Matchers
    with MockSugar
    with PrivateMethodTester {

  val startStatusEvent = StatusEvent
    .newBuilder()
    .setType("STATUS_CHANGED")
    .setEventTime(Timestamp.newBuilder().setSeconds(1).build())
    .setDescription("Job state is set from SCHEDULED to RUNNING for job...")
    .build()

  val endStatusEvent = StatusEvent
    .newBuilder()
    .setType("STATUS_CHANGED")
    .setEventTime(Timestamp.newBuilder().setSeconds(2).build())
    .setDescription("Job state is set from RUNNING to SOME_OTHER_STATUS for job...")
    .build()

  val schedulingStatusEvent = StatusEvent
    .newBuilder()
    .setType("STATUS_CHANGED")
    .setEventTime(Timestamp.newBuilder().setSeconds(3).build())
    .setDescription("Job state is set from QUEUED to SCHEDULED for job...")
    .build()

  val preemptionError = "Job state is set from SCHEDULED to FAILED for job projects/....Job failed due to task " +
    "failure. Specifically, task with index 0 failed due to the following task event: \"Task state is updated " +
    "from PENDING to FAILED on zones/... due to Spot VM preemption with exit code 50001.\""
  val preemptionStatusEvent = StatusEvent
    .newBuilder()
    .setType("STATUS_CHANGED")
    .setEventTime(Timestamp.newBuilder().setSeconds(4).build())
    .setDescription(preemptionError)
    .build()

  def setupBatchClient(machineType: String = "n1-standard-1",
                       location: String = "regions/us-central1",
                       jobState: State = JobStatus.State.SUCCEEDED,
                       events: List[StatusEvent]
  ): BatchServiceClient = {
    val instancePolicy = InstancePolicy
      .newBuilder()
      .setMachineType(machineType)
      .setProvisioningModel(ProvisioningModel.PREEMPTIBLE)
      .build()

    val allocationPolicy = AllocationPolicy
      .newBuilder()
      .setLocation(LocationPolicy.newBuilder().addAllowedLocations(location))
      .addInstances(InstancePolicyOrTemplate.newBuilder().setPolicy(instancePolicy))
      .build()

    val jobStatus = JobStatus
      .newBuilder()
      .setState(jobState)
      .addAllStatusEvents(events.asJava)
      .build()

    val job = Job.newBuilder().setAllocationPolicy(allocationPolicy).setStatus(jobStatus).build()

    val mockClient = mock[BatchServiceClient]
    doReturn(job).when(mockClient).getJob(any[GetJobRequest])
    doReturn(job).when(mockClient).getJob(any[String])

    mockClient
  }

  behavior of "BatchRequestExecutor"

  it should "create a schedule status event correctly" in {

    val mockClient = setupBatchClient(jobState = JobStatus.State.QUEUED, events = List(schedulingStatusEvent))
    val batchRequestExecutor = new BatchRequestExecutor.CloudImpl(BatchServiceSettings.newBuilder().build())

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val internalGetHandler = PrivateMethod[BatchApiResponse.StatusQueried](Symbol("internalGetHandler"))
    val result = batchRequestExecutor invokePrivate internalGetHandler(mockClient, GetJobRequest.newBuilder().build())

    // Verify the event
    result.status match {
      case RunStatus.Initializing(events, _) =>
        events.length shouldBe 1
        events.map(_.name).head shouldBe "Job state is set from QUEUED to SCHEDULED for job..."
        events.map(_.offsetDateTime.toString).head shouldBe "1970-01-01T00:00:03Z"
      case _ => fail("Expected RunStatus.Initializing with events")
    }
  }

  it should "handle an immediate preemption correctly" in {

    val mockClient = setupBatchClient(jobState = JobStatus.State.FAILED, events = List(preemptionStatusEvent))
    val batchRequestExecutor = new BatchRequestExecutor.CloudImpl(BatchServiceSettings.newBuilder().build())

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val internalGetHandler = PrivateMethod[BatchApiResponse.StatusQueried](Symbol("internalGetHandler"))
    val result = batchRequestExecutor invokePrivate internalGetHandler(mockClient, GetJobRequest.newBuilder().build())

    // Verify the event
    result.status match {
      case RunStatus.Preempted(_, events, _) =>
        events.length shouldBe 1
        events.map(_.name).head shouldBe preemptionError
        events.map(_.offsetDateTime.toString).head shouldBe "1970-01-01T00:00:04Z"
      case _ => fail("Expected RunStatus.Preempted with events")
    }
  }

  it should "create instantiatedVmInfo correctly" in {

    val mockClient =
      setupBatchClient(jobState = JobStatus.State.RUNNING, events = List(startStatusEvent, endStatusEvent))
    // Create the BatchRequestExecutor
    val batchRequestExecutor = new BatchRequestExecutor.CloudImpl(BatchServiceSettings.newBuilder().build())

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val internalGetHandler = PrivateMethod[BatchApiResponse.StatusQueried](Symbol("internalGetHandler"))
    val result = batchRequestExecutor invokePrivate internalGetHandler(mockClient, GetJobRequest.newBuilder().build())

    // Verify the instantiatedVmInfo
    result.status match {
      case RunStatus.Running(_, Some(instantiatedVmInfo)) =>
        instantiatedVmInfo.region shouldBe "us-central1"
        instantiatedVmInfo.machineType shouldBe "n1-standard-1"
        instantiatedVmInfo.preemptible shouldBe true
      case _ => fail("Expected RunStatus.Running with instantiatedVmInfo")
    }
  }

  it should "create instantiatedVmInfo correctly with different location info" in {

    val mockClient = setupBatchClient(location = "zones/us-central1-a",
                                      jobState = JobStatus.State.RUNNING,
                                      events = List(startStatusEvent, endStatusEvent)
    )

    // Create the BatchRequestExecutor
    val batchRequestExecutor = new BatchRequestExecutor.CloudImpl(BatchServiceSettings.newBuilder().build())

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val internalGetHandler = PrivateMethod[BatchApiResponse.StatusQueried](Symbol("internalGetHandler"))
    val result = batchRequestExecutor invokePrivate internalGetHandler(mockClient, GetJobRequest.newBuilder().build())

    // Verify the instantiatedVmInfo
    result.status match {
      case RunStatus.Running(_, Some(instantiatedVmInfo)) =>
        instantiatedVmInfo.region shouldBe "us-central1-a"
        instantiatedVmInfo.machineType shouldBe "n1-standard-1"
        instantiatedVmInfo.preemptible shouldBe true
      case _ => fail("Expected RunStatus.Running with instantiatedVmInfo")
    }
  }

  it should "create instantiatedVmInfo correctly with missing location info" in {

    val mockClient =
      setupBatchClient(jobState = JobStatus.State.RUNNING, events = List(startStatusEvent, endStatusEvent))

    // Create the BatchRequestExecutor
    val batchRequestExecutor = new BatchRequestExecutor.CloudImpl(BatchServiceSettings.newBuilder().build())

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val internalGetHandler = PrivateMethod[BatchApiResponse.StatusQueried](Symbol("internalGetHandler"))
    val result = batchRequestExecutor invokePrivate internalGetHandler(mockClient, GetJobRequest.newBuilder().build())

    // Verify the instantiatedVmInfo
    result.status match {
      case RunStatus.Running(_, Some(instantiatedVmInfo)) =>
        instantiatedVmInfo.region shouldBe "us-central1"
        instantiatedVmInfo.machineType shouldBe "n1-standard-1"
        instantiatedVmInfo.preemptible shouldBe true
      case _ => fail("Expected RunStatus.Running with instantiatedVmInfo")
    }
  }

  it should "send vmStartTime and vmEndTime metadata info when a workflow succeeds" in {

    val mockClient = setupBatchClient(events = List(startStatusEvent, endStatusEvent))

    // Create the BatchRequestExecutor
    val batchRequestExecutor = new BatchRequestExecutor.CloudImpl(BatchServiceSettings.newBuilder().build())

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val internalGetHandler = PrivateMethod[BatchApiResponse.StatusQueried](Symbol("internalGetHandler"))
    val result = batchRequestExecutor invokePrivate internalGetHandler(mockClient, GetJobRequest.newBuilder().build())

    // Verify the events
    result.status match {
      case RunStatus.Success(events, _) =>
        val eventNames = events.map(_.name)
        val eventTimes = events.map(_.offsetDateTime.toString)
        eventNames should contain allOf ("vmStartTime", "vmEndTime")
        eventTimes should contain allOf ("1970-01-01T00:00:01Z", "1970-01-01T00:00:02Z")
      case _ => fail("Expected RunStatus.Success with events")
    }
  }

  it should "send vmStartTime and vmEndTime metadata info along with other events when a workflow fails" in {
    val mockClient =
      setupBatchClient(jobState = JobStatus.State.FAILED, events = List(startStatusEvent, endStatusEvent))

    // Create the BatchRequestExecutor
    val batchRequestExecutor = new BatchRequestExecutor.CloudImpl(BatchServiceSettings.newBuilder().build())

    // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
    val internalGetHandler = PrivateMethod[BatchApiResponse.StatusQueried](Symbol("internalGetHandler"))
    val result = batchRequestExecutor invokePrivate internalGetHandler(mockClient, GetJobRequest.newBuilder().build())

    // Verify the events
    result.status match {
      case RunStatus.Failed(_, events, _) =>
        val eventNames = events.map(_.name)
        val eventTimes = events.map(_.offsetDateTime.toString)
        println(eventNames)
        eventNames should contain allOf ("vmStartTime", "vmEndTime")
        eventNames should contain allOf ("Job state is set from RUNNING to SOME_OTHER_STATUS for job...", "Job state is set from SCHEDULED to RUNNING for job...")
        eventTimes should contain allOf ("1970-01-01T00:00:01Z", "1970-01-01T00:00:02Z")
      case _ => fail("Expected RunStatus.Success with events")
    }
  }

}
