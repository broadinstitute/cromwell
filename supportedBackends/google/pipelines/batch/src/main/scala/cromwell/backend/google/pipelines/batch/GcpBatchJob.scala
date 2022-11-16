package cromwell.backend.google.pipelines.batch

import akka.actor.ActorRef
//import akka.actor.{Actor, ActorRef}
import com.google.cloud.batch.v1.AllocationPolicy.{InstancePolicy, InstancePolicyOrTemplate}
import com.google.cloud.batch.v1.LogsPolicy.Destination
import com.google.cloud.batch.v1.{AllocationPolicy, BatchServiceClient, ComputeResource, CreateJobRequest, Job, LogsPolicy, Runnable, TaskGroup, TaskSpec}
import com.google.protobuf.Duration
import com.google.cloud.batch.v1.Runnable.Container
import cromwell.backend.BackendJobDescriptor

import java.util.concurrent.TimeUnit
//import scala.concurrent.Promise



final case class GcpBatchJob(jobDescriptor: BackendJobDescriptor)  {

  val projectId = "batch-testing-350715"
  val region = "us-central1"
  val jobName = "hello-dan-container-1"


  val batchServiceClient = BatchServiceClient
    .create


  def callClient: Unit = {

    val runnable = Runnable
      .newBuilder
      .setContainer((Container
        .newBuilder
        .setImageUri("gcr.io/google-containers/busybox")
        .setEntrypoint("/bin/sh")
        .addCommands("-c")
        .addCommands("echo Hello World!")
        .build
        ))
      .build

    val computeResource = ComputeResource
      .newBuilder
      .setCpuMilli(2000)
      .setMemoryMib(200)
      .build

    //define task spec
    val task = TaskSpec
      .newBuilder
      .addRunnables(runnable)
      .setComputeResource(computeResource)
      .setMaxRetryCount(2)
      .setMaxRunDuration(Duration
        .newBuilder
        .setSeconds(3600)
        .build) // Jobs can be divided into tasks. In this case, we have only one task.

    // Tasks are grouped inside a job using TaskGroups.
    val taskGroup = TaskGroup
      .newBuilder
      .setTaskCount(4)
      .setTaskSpec(task)
      .build

    // Instance Policy which also defines machine type
    val instancePolicy = InstancePolicy
      .newBuilder
      .setMachineType("e2-standard-4")
      .build

    // Allocation policy
    val allocationPolicy = AllocationPolicy
      .newBuilder
      .addInstances(InstancePolicyOrTemplate
        .newBuilder
        .setPolicy(instancePolicy)
        .build)
      .build

    // Job.  Includes definitions from earlier
    val job = Job
      .newBuilder
      .addTaskGroups(taskGroup)
      .setAllocationPolicy(allocationPolicy)
      .putLabels("env", "testing")
      .putLabels("type", "script")
      .setLogsPolicy(LogsPolicy
        .newBuilder
        .setDestination(Destination
          .CLOUD_LOGGING)
        .build)

    //set parent
    val parent = (String
      .format("projects/%s/locations/%s", projectId, region))

    //Create job request.  Run
    val createJobRequest = CreateJobRequest
      .newBuilder
      .setParent(parent)
      .setJob(job)
      .setJobId(jobName)
      .build()

    val result = batchServiceClient.createJobCallable.futureCall(createJobRequest).get(3, TimeUnit .MINUTES)

    result.getName


  }

}


