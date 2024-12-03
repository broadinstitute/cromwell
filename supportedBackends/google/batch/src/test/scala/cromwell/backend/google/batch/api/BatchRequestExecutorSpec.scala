package cromwell.backend.google.batch.api.request

import akka.actor.ActorSystem
import akka.testkit.TestKit
import com.google.cloud.batch.v1.{AllocationPolicy, BatchServiceClient, BatchServiceSettings, GetJobRequest, Job, JobStatus}
import com.google.cloud.batch.v1.AllocationPolicy.{InstancePolicy, InstancePolicyOrTemplate, LocationPolicy, ProvisioningModel}
import common.mock.MockSugar
import cromwell.backend.google.batch.api.BatchApiResponse
import cromwell.backend.google.batch.models.RunStatus
import org.mockito.ArgumentMatchers.any
import org.mockito.Mockito.doReturn
import org.scalatest.PrivateMethodTester
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class BatchRequestExecutorSpec
    extends TestKit(ActorSystem("BatchRequestExecutorSpec"))
    with AnyFlatSpecLike
    with Matchers
    with MockSugar
    with PrivateMethodTester {

  behavior of "BatchRequestExecutor"

  it should "create instantiatedVmInfo correctly" in {

    val instancePolicy = InstancePolicy
      .newBuilder()
      .setMachineType("n1-standard-1")
      .setProvisioningModel(ProvisioningModel.PREEMPTIBLE)
      .build()

    val allocationPolicy = AllocationPolicy
      .newBuilder()
      .setLocation(LocationPolicy.newBuilder().addAllowedLocations("regions/us-central1"))
      .addInstances(InstancePolicyOrTemplate.newBuilder().setPolicy(instancePolicy))
      .build()

    val jobStatus = JobStatus.newBuilder().setState(JobStatus.State.RUNNING).build()

    val job = Job.newBuilder().setAllocationPolicy(allocationPolicy).setStatus(jobStatus).build()

    val mockClient = mock[BatchServiceClient]
    doReturn(job).when(mockClient).getJob(any[GetJobRequest])
    doReturn(job).when(mockClient).getJob(any[String])

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
}
