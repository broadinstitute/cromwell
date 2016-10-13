package cromwell.backend.impl.aws

import akka.actor.Props
import com.amazonaws.auth.AWSCredentials
import com.amazonaws.services.ecs.AmazonECSAsyncClient
import cromwell.backend.BackendJobExecutionActor.{BackendJobExecutionResponse, SucceededResponse}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobExecutionActor}
import com.amazonaws.services.ecs.model._
import cromwell.backend.impl.aws.util.AwsSdkAsyncHandler
import cromwell.backend.impl.aws.util.AwsSdkAsyncHandler.AwsSdkAsyncResult
import cromwell.backend.wdl.OnlyPureFunctions
import net.ceedubs.ficus.Ficus._

import scala.collection.JavaConverters._
import scala.concurrent.duration.Duration
import scala.concurrent.{Await, Future}

class AwsJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                           override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {

  val awsAccessKeyId = configurationDescriptor.backendConfig.as[String]("accessKeyId")
  val awsSecretKey = configurationDescriptor.backendConfig.as[String]("secretKey")

  val clusterName = "ecs-t2micro-cluster"

  val credentials = new AWSCredentials {
    override def getAWSAccessKeyId: String = awsAccessKeyId
    override def getAWSSecretKey: String = awsSecretKey
  }
  val ecsAsyncClient = new AmazonECSAsyncClient(credentials)

  override def execute: Future[BackendJobExecutionResponse] = {

    val commandOverride = new ContainerOverride().withName("simple-app").withCommand(jobDescriptor.call.instantiateCommandLine(Map.empty, OnlyPureFunctions, identity).get)

    val runRequest: RunTaskRequest = new RunTaskRequest()
      .withCluster(clusterName)
      .withCount(1)
      .withTaskDefinition("ubuntuTask:1")
      .withOverrides(new TaskOverride().withContainerOverrides(commandOverride))

    val submitResultHandler = new AwsSdkAsyncHandler[RunTaskRequest, RunTaskResult]()
    val _ = ecsAsyncClient.runTaskAsync(runRequest, submitResultHandler)

    submitResultHandler.future map {
      case AwsSdkAsyncResult(_, result) =>
        log.info("AWS submission completed:\n{}", result.toString)
        val taskArn= result.getTasks.asScala.head.getTaskArn
        val taskDescription = waitUntilDone(taskArn)

        log.info("AWS task completed!\n{}", taskDescription.toString)
        SucceededResponse(jobDescriptor.key, Option(0), Map.empty, None, Seq.empty)
    }
  }

  private def waitUntilDone(taskArn: String): Task = {
    val describeTasksRequest = new DescribeTasksRequest().withCluster(clusterName).withTasks(List(taskArn).asJava)

    val resultHandler = new AwsSdkAsyncHandler[DescribeTasksRequest, DescribeTasksResult]()
    val _ = ecsAsyncClient.describeTasksAsync(describeTasksRequest, resultHandler)

    val describedTasks = Await.result(resultHandler.future, Duration.Inf)
    val taskDescription = describedTasks.result.getTasks.asScala.head
    if (taskDescription.getLastStatus == DesiredStatus.STOPPED.toString) {
      taskDescription
    } else {
      log.info(s"Still waiting for completion. Last known status: {}", taskDescription.getLastStatus)
      Thread.sleep(2000)
      waitUntilDone(taskArn)
    }
  }
}

object AwsJobExecutionActor {
  def props(jobDescriptor: BackendJobDescriptor,
            configurationDescriptor: BackendConfigurationDescriptor): Props = Props(new AwsJobExecutionActor(jobDescriptor, configurationDescriptor))
}
