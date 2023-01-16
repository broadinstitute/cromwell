package cromwell.backend.google.pipelines.batch
import com.google.cloud.batch.v1.AllocationPolicy.{InstancePolicy, InstancePolicyOrTemplate, LocationPolicy}
import com.google.cloud.batch.v1.Runnable.Container
import com.google.cloud.batch.v1.{AllocationPolicy, BatchServiceClient, ComputeResource, CreateJobRequest, Job, LogsPolicy, Runnable, TaskGroup, TaskSpec}
import cromwell.backend.google.pipelines.batch.GcpBatchBackendSingletonActor.BatchRequest
import com.google.protobuf.Duration
import com.google.cloud.batch.v1.LogsPolicy.Destination

import java.util.concurrent.TimeUnit

final case class GcpBatchJob (
                             jobSubmission: BatchRequest,
                             cpu: Long,
                             memory: Long,
                             machineType: String,
                             dockerImage: String,
                            ) {

  // VALUES HERE
  val entryPoint = "/bin/sh"
  val retryCount = 2
  val durationInSeconds: Long = 3600
  val taskCount: Long = 1

  lazy val batchServiceClient = BatchServiceClient.create

  lazy val parent = (String
    .format("projects/%s/locations/%s", jobSubmission
      .projectId, jobSubmission
      .region))
  private def createRunnable(dockerImage: String, entryPoint: String): Runnable = {
    val runnable = Runnable.newBuilder.setContainer((Container.newBuilder.setImageUri(dockerImage).setEntrypoint(entryPoint).addCommands("-c").addCommands("echo Hello World!").build)).build
    runnable
  }

  private def createInstancePolicy = {
    val instancePolicy = InstancePolicy
      .newBuilder
      .setMachineType(machineType)
      .build
    instancePolicy
  }

  def submitJob(): Unit = {

    try {
      //val runnable = Runnable.newBuilder.setContainer((Container.newBuilder.setImageUri(dockerImage).setEntrypoint(entryPoint).addCommands("-c").addCommands("echo Hello World!").build)).build
      val runnable = createRunnable(dockerImage = dockerImage, entryPoint = entryPoint)
      val computeResource = ComputeResource
        .newBuilder
        .setCpuMilli(cpu)
        .setMemoryMib(memory)
        .build
      val task = TaskSpec
        .newBuilder
        .addRunnables(runnable)
        .setComputeResource(computeResource)
        .setMaxRetryCount(retryCount)
        .setMaxRunDuration(Duration
          .newBuilder
          .setSeconds(durationInSeconds)
          .build)
      val taskGroup = TaskGroup
        .newBuilder
        .setTaskCount(taskCount)
        .setTaskSpec(task)
        .build
      val instancePolicy = createInstancePolicy
      val locationPolicy = LocationPolicy.newBuilder.addAllowedLocations("regions/us-south1")
      //val locationPolicy = LocationPolicy.newBuilder.setAllowedLocations(0, "regions/us-south-1")
      val allocationPolicy = AllocationPolicy
        .newBuilder
        .setLocation(locationPolicy)
        .addInstances(InstancePolicyOrTemplate
          .newBuilder
          .setPolicy(instancePolicy)
          .build)
        .build
      val job = Job
        .newBuilder
        .addTaskGroups(taskGroup)
        .setAllocationPolicy(allocationPolicy)
        .putLabels("submitter", "cromwell") // label to signify job submitted by cromwell for larger tracking purposes within GCP batch
        .putLabels("env", "testing")
        .putLabels("type", "script")
        .setLogsPolicy(LogsPolicy
          .newBuilder
          .setDestination(Destination
            .CLOUD_LOGGING)
          .build)

      val createJobRequest = CreateJobRequest
        .newBuilder
        .setParent(parent)
        .setJob(job)
        .setJobId(jobSubmission
          .jobName)
        .build()
      batchServiceClient
        .createJobCallable
        .futureCall(createJobRequest)
        .get(3, TimeUnit
          .MINUTES)
      println("job submitted")

    }
    catch  {
      case _: Throwable => println("Job failed")
    }

  }
}