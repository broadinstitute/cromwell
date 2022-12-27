package cromwell.backend.google.pipelines.batch

import akka.actor.{Actor, ActorLogging, Props}

//import scala.concurrent.duration.DurationInt
//import akka.actor.{Actor, ActorSystem, ActorLogging, Props}
//import akka.actor.{Actor, ActorSystem, ActorLogging, ActorRef, Props}
import com.google.cloud.batch.v1.AllocationPolicy.{InstancePolicy, InstancePolicyOrTemplate}
import com.google.cloud.batch.v1.LogsPolicy.Destination
import com.google.cloud.batch.v1.Runnable.Container
import com.google.cloud.batch.v1.{AllocationPolicy, BatchServiceClient, ComputeResource, CreateJobRequest, Job, LogsPolicy, Runnable, TaskGroup, TaskSpec}
import com.google.protobuf.Duration
import scala.concurrent.ExecutionContext
//import cromwell.core.Dispatcher.BackendDispatcher
//import cats.effect.{IO, Timer}

import java.util.concurrent.TimeUnit

object GcpBatchBackendSingletonActor {
  def props(name: String) = Props(new GcpBatchBackendSingletonActor(name))

  case class BatchRequest(projectId: String, region: String, jobName: String)
}

class GcpBatchBackendSingletonActor (name: String) extends Actor with ActorLogging {

  import GcpBatchBackendSingletonActor._

  implicit val ec: ExecutionContext = context.dispatcher

  def receive: Receive = {

    case jobSubmission: BatchRequest =>
      val batchServiceClient = BatchServiceClient.create
      val runnable = Runnable.newBuilder.setContainer((Container.newBuilder.setImageUri("gcr.io/google-containers/busybox").setEntrypoint("/bin/sh").addCommands("-c").addCommands("echo Hello World!").build)).build
      val computeResource = ComputeResource.newBuilder.setCpuMilli(2000).setMemoryMib(200).build
      val task = TaskSpec.newBuilder.addRunnables(runnable).setComputeResource(computeResource).setMaxRetryCount(2).setMaxRunDuration(Duration.newBuilder.setSeconds(3600).build)
      val taskGroup = TaskGroup.newBuilder.setTaskCount(1).setTaskSpec(task).build
      val instancePolicy = InstancePolicy.newBuilder.setMachineType("e2-standard-4").build
      val allocationPolicy = AllocationPolicy.newBuilder.addInstances(InstancePolicyOrTemplate.newBuilder.setPolicy(instancePolicy).build).build
      val job = Job.newBuilder.addTaskGroups(taskGroup).setAllocationPolicy(allocationPolicy).putLabels("env", "testing").putLabels("type", "script").setLogsPolicy(LogsPolicy.newBuilder.setDestination(Destination.CLOUD_LOGGING).build)
      val parent = (String.format("projects/%s/locations/%s", jobSubmission.projectId, jobSubmission.region))
      val createJobRequest = CreateJobRequest.newBuilder.setParent(parent).setJob(job).setJobId(jobSubmission.jobName).build()
      val result = batchServiceClient.createJobCallable.futureCall(createJobRequest).get(3, TimeUnit.MINUTES)
      println(result.getName)

    case other =>
      log.error("Unknown message to GCP Batch Singleton Actor: {}. Dropping it.", other)

  }

}


