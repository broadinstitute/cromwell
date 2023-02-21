package cromwell.backend.google.pipelines.batch
import com.google.api.gax.rpc.{FixedHeaderProvider, HeaderProvider}
import com.google.cloud.batch.v1.{AllocationPolicy, BatchServiceClient, BatchServiceSettings, ComputeResource, CreateJobRequest, Job, LogsPolicy, Runnable, TaskGroup, TaskSpec}
//import com.google.cloud.batch.v1.{AllocationPolicy, BatchServiceClient, BatchServiceSettings, ComputeResource, CreateJobRequest, GetJobRequest, Job, JobName, LogsPolicy, Runnable, TaskGroup, TaskSpec}
import com.google.cloud.batch.v1.AllocationPolicy.{InstancePolicy, InstancePolicyOrTemplate, LocationPolicy}
import com.google.cloud.batch.v1.Runnable.Container
import cromwell.backend.google.pipelines.batch.GcpBatchBackendSingletonActor.BatchRequest
import com.google.protobuf.Duration
import com.google.cloud.batch.v1.LogsPolicy.Destination
import com.google.common.collect.ImmutableMap

import java.util.concurrent.TimeUnit
import scala.util.Try

final case class GcpBatchJob (
                             jobSubmission: BatchRequest,
                             cpu: Long,
                             //cpuPlatform: String,
                             memory: Long,
                             machineType: String,
                             runtimeAttributes: GcpBatchRuntimeAttributes
                            ) {

  // VALUES HERE
  val entryPoint = "/bin/sh"
  val retryCount = 2
  val durationInSeconds: Long = 3600
  val taskCount: Long = 1

  // set user agent
  val user_agent_header = "user-agent"
  val customUserAgentValue = "cromwell"
  lazy val headerProvider: HeaderProvider = FixedHeaderProvider
    .create(ImmutableMap
      .of(user_agent_header, customUserAgentValue))

  lazy val batchSettings = BatchServiceSettings.newBuilder.setHeaderProvider(headerProvider).build

  lazy val batchServiceClient = BatchServiceClient.create(batchSettings)

  lazy val parent = (String
    .format("projects/%s/locations/%s", jobSubmission
      .projectId, jobSubmission
      .region))

  private val cpuPlatform =  runtimeAttributes.cpuPlatform.getOrElse("")

  println(cpuPlatform)
  private def createRunnable(dockerImage: String, entryPoint: String): Runnable = {
    val runnable = Runnable.newBuilder.setContainer((Container.newBuilder.setImageUri(dockerImage).setEntrypoint(entryPoint).addCommands("-c").addCommands("echo Hello World!").build)).build
    runnable
  }

  private def createInstancePolicy = {
    val instancePolicy = InstancePolicy
      .newBuilder
      .setMachineType(machineType)
      .setMinCpuPlatform(cpuPlatform)
      .build
    instancePolicy
  }

  def submitJob(): Unit = {

    try {
      val runnable = createRunnable(dockerImage = runtimeAttributes.dockerImage, entryPoint = entryPoint)
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
        .putLabels("cromwell-workflow-id", jobSubmission.workflowId.toString) // label to make it easier to match Cromwell workflows with multiple GCP batch jobs
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
      val result = batchServiceClient
        .createJobCallable
        .futureCall(createJobRequest)
        .get(5, TimeUnit
          .SECONDS)
      println("job submitted")

      batchServiceClient.close()

      println(result.getName)


    }
    catch  {
      case _: Throwable => println("Job failed")
    }

  }


  def jobGetRequest(jobId: String) = {
    val gcpBatchPoll = new GcpBatchJobGetRequest
    gcpBatchPoll.GetJob(jobId)

  }
  def status(jobId: String): Try[RunStatus] = for {
    _ <- Try(jobGetRequest(jobId).toString)
    //runStatus <- RunStatus.fromJobStatus(jobId)
    runStatus <- RunStatus.testJobStatus(jobId)
  } yield runStatus

}