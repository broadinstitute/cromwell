package cromwell.backend.google.pipelines.batch
import com.google.cloud.batch.v1.AllocationPolicy.{InstancePolicy, InstancePolicyOrTemplate}
import com.google.cloud.batch.v1.Runnable.Container
import com.google.cloud.batch.v1.{AllocationPolicy, BatchServiceClient, ComputeResource, CreateJobRequest, Job, LogsPolicy, Runnable, TaskGroup, TaskSpec}
import cromwell.backend.google.pipelines.batch.GcpBatchBackendSingletonActor.BatchRequest
import com.google.protobuf.Duration
import com.google.cloud.batch.v1.LogsPolicy.Destination
import java.util.concurrent.TimeUnit

final case class GcpBatchJob(
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
  val parentPath = "projects/%s/locations/%s"

  def submitJob(): Unit = {
    try {
      val batchServiceClient = BatchServiceClient.create
      val runnable = Runnable.newBuilder.setContainer((Container.newBuilder.setImageUri(dockerImage).setEntrypoint(entryPoint).addCommands("-c").addCommands("echo Hello World!").build)).build
      val computeResource = ComputeResource.newBuilder.setCpuMilli(cpu).setMemoryMib(memory).build
      val task = TaskSpec.newBuilder.addRunnables(runnable).setComputeResource(computeResource).setMaxRetryCount(retryCount).setMaxRunDuration(Duration.newBuilder.setSeconds(durationInSeconds).build)
      val taskGroup = TaskGroup.newBuilder.setTaskCount(taskCount).setTaskSpec(task).build
      val instancePolicy = InstancePolicy.newBuilder.setMachineType(machineType).build
      val allocationPolicy = AllocationPolicy.newBuilder.addInstances(InstancePolicyOrTemplate.newBuilder.setPolicy(instancePolicy).build).build
      val job = Job.newBuilder.addTaskGroups(taskGroup).setAllocationPolicy(allocationPolicy).putLabels("env", "testing").putLabels("type", "script").setLogsPolicy(LogsPolicy.newBuilder.setDestination(Destination.CLOUD_LOGGING).build)
      val parent = (String.format(parentPath, jobSubmission.projectId, jobSubmission.region))
      val createJobRequest = CreateJobRequest.newBuilder.setParent(parent).setJob(job).setJobId(jobSubmission.jobName).build()
      val result = batchServiceClient.createJobCallable.futureCall(createJobRequest).get(3, TimeUnit.MINUTES)
    } catch {
      case _: Throwable => println("OOPS") // clean up
    }

  }
}