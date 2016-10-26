package cromwell.backend.impl.aws

import java.nio.file.Path

import akka.actor.{Actor, ActorLogging, ActorRef}
import better.files.File
import com.amazonaws.auth.{AWSCredentials, BasicAWSCredentials}
import com.amazonaws.services.ecs.model._
import com.amazonaws.services.ecs.{AmazonECSAsync, AmazonECSAsyncClient}
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.async.AsyncBackendJobExecutionActor.{ExecutionMode, Recover}
import cromwell.backend.async.{AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, PendingExecutionHandle, SuccessfulExecutionHandle}
import cromwell.backend.validation.{ContinueOnReturnCode, ContinueOnReturnCodeFlag}
import cromwell.backend.wdl.{Command, OnlyPureFunctions}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor, BackendJobLifecycleActor}
import cromwell.core.retry.SimpleExponentialBackoff
import net.ceedubs.ficus.Ficus._
import wdl4s.expression.WdlFunctions
import wdl4s.values.WdlValue
import wdl4s.{Task => _, _}

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Failure, Success, Try}

class AwsAsyncJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                override val completionPromise: Promise[BackendJobExecutionResponse],
                                val configurationDescriptor: BackendConfigurationDescriptor,
                                val serviceRegistryActor: ActorRef)
  extends Actor with ActorLogging with BackendJobLifecycleActor with AsyncBackendJobExecutionActor {

  // Copypasta

  override def retryable = false

  override def executeOrRecoverBackOff = SimpleExponentialBackoff(
    initialInterval = 3.seconds, maxInterval = 20.seconds, multiplier = 1.1)

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    Future.fromTry(Try {
      mode match {
        case Recover(jobId: BackendJobId) => recover(jobId)
        case _ => execute()
      }
    })
  }

  override lazy val pollBackOff = SimpleExponentialBackoff(
    initialInterval = 30.seconds, maxInterval = 600.seconds, multiplier = 1.1)

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    Future.fromTry(Try {
      previous match {
        case handle: PendingExecutionHandle[
          BackendJobId@unchecked, BackendRunInfo@unchecked, BackendRunStatus@unchecked] =>

          jobLogger.debug(s"$tag Polling Job ${handle.pendingJob}")
          Try(pollStatus(handle.pendingJob)) match {
            case Success(backendRunStatus) => updateExecutionHandleSuccess(handle, backendRunStatus)
            case Failure(throwable) => updateExecutionHandleFailure(handle, throwable)
          }
        case f: FailedNonRetryableExecutionHandle => f
        case s: SuccessfulExecutionHandle => s
        case badHandle => throw new IllegalArgumentException(s"Unexpected execution handle: $badHandle")
      }
    })
  }

  private def updateExecutionHandleSuccess(oldHandle: BackendPendingExecutionHandle,
                                           status: BackendRunStatus): ExecutionHandle = {
    val previousStatus = oldHandle.previousStatus
    if (!(previousStatus contains status)) {
      // If this is the first time checking the status, we log the transition as '-' to 'currentStatus'. Otherwise
      // just use the state names.
      val prevStateName = previousStatus.map(_.toString).getOrElse("-")
      jobLogger.info(s"$tag Status change from $prevStateName to $status")
      tellMetadata(Map("backendStatus" -> status))
    }

    status match {
      case _ if isTerminal(status) =>
        val metadata = getTerminalMetadata(status)
        tellMetadata(metadata)
        executionResult(status, oldHandle)
      case s => oldHandle.copy(previousStatus = Option(s)) // Copy the current handle with updated previous status.
    }
  }

  private def updateExecutionHandleFailure(oldHandle: BackendPendingExecutionHandle,
                                           throwable: Throwable): ExecutionHandle = {
    throwable match {
      case exception: Exception =>
        val handler: PartialFunction[Exception, ExecutionHandle] =
          customPollStatusFailure(oldHandle) orElse {
            case exception: Exception =>
              // Log exceptions and return the original handle to try again.
              jobLogger.warn(s"Caught exception, retrying", exception)
              oldHandle
          }
        handler(exception)
      case error: Error => throw error // JVM-ending calamity.
      case _: Throwable =>
        // Someone has subclassed or instantiated Throwable directly. Kill the job. They should be using an Exception.
        FailedNonRetryableExecutionHandle(throwable)
    }
  }

  private def tellMetadata(metadataKeyValues: Map[String, Any]): Unit = {
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(jobDescriptor.workflowDescriptor.id, Option(jobDescriptor.key), metadataKeyValues)
  }

  private def executionResult(status: BackendRunStatus, handle: BackendPendingExecutionHandle)
                             (implicit ec: ExecutionContext): ExecutionHandle = {
    try {
      if (isSuccess(status)) {

        lazy val stderrLength: Long = File(remoteStrErrPath).size
        lazy val returnCode: Try[Int] = returnCodeContents.map(_.trim.toInt)
        status match {
          case _ if failOnStderr && stderrLength.intValue > 0 =>
            // returnCode will be None if it couldn't be downloaded/parsed, which will yield a null in the DB
            FailedNonRetryableExecutionHandle(new RuntimeException(
              s"execution failed: stderr has length $stderrLength"), returnCode.toOption)
          case _ if returnCodeContents.isFailure =>
            val exception = returnCode.failed.get
            jobLogger.warn(s"could not download return code file, retrying", exception)
            // Return handle to try again.
            handle
          case _ if returnCode.isFailure =>
            FailedNonRetryableExecutionHandle(new RuntimeException(
              s"execution failed: could not parse return code as integer: ${returnCodeContents.get}"))
          case _ if !continueOnReturnCode.continueFor(returnCode.get) =>
            val badReturnCodeMessage = s"Call ${jobDescriptor.key}: return code was ${returnCode.getOrElse("(none)")}"
            FailedNonRetryableExecutionHandle(new RuntimeException(badReturnCodeMessage), returnCode.toOption)
          case _ =>
            handleSuccess(status, handle, returnCode.get)
        }

      } else {
        handleFailure(status, handle)
      }
    } catch {
      case e: Exception =>
        jobLogger.warn("Caught exception processing job result, retrying", e)
        // Return the original handle to try again.
        handle
    }
  }

  private lazy val instantiatedCommand = Command.instantiate(
    jobDescriptor, commandLineFunctions, commandLinePreProcessor, commandLineValueMapper).get

  private lazy val tag = s"${this.getClass.getSimpleName} [UUID(${workflowId.shortString}):${jobDescriptor.key.tag}]"

  // aws backend specific

  type BackendJobId = AwsJobId
  type BackendRunInfo = Any
  type BackendRunStatus = AwsRunStatus
  type BackendPendingExecutionHandle = PendingExecutionHandle[BackendJobId, BackendRunInfo, BackendRunStatus]

  private def commandLineFunctions: WdlFunctions[WdlValue] = OnlyPureFunctions

  private def commandLinePreProcessor: EvaluatedTaskInputs => Try[EvaluatedTaskInputs] = Success.apply

  private def commandLineValueMapper: WdlValue => WdlValue = identity

  private def remoteStrErrPath: Path = ???

  private lazy val returnCodeContents: Try[String] = Try("0")

  private def continueOnReturnCode: ContinueOnReturnCode = ContinueOnReturnCodeFlag(false)

  private def failOnStderr: Boolean = false

  private val awsAccessKeyId: String = configurationDescriptor.backendConfig.as[String]("accessKeyId")
  private val awsSecretKey: String = configurationDescriptor.backendConfig.as[String]("secretKey")
  private val clusterName: String =
    configurationDescriptor.backendConfig.getOrElse("clusterName", "ecs-t2micro-cluster")

  private val credentials: AWSCredentials = new BasicAWSCredentials(awsAccessKeyId, awsSecretKey)
  private val ecsAsyncClient: AmazonECSAsync = new AmazonECSAsyncClient(credentials)

  private def execute(): ExecutionHandle = {
    val commandOverride = new ContainerOverride()
      .withName("simple-app")
      .withCommand(instantiatedCommand)

    val taskOverride = new TaskOverride()
      .withContainerOverrides(commandOverride)

    val runTaskRequest = new RunTaskRequest()
      .withCluster(clusterName)
      .withCount(1)
      .withTaskDefinition("ubuntuTask:1")
      .withOverrides(taskOverride)

    val runTaskResult = ecsAsyncClient.runTask(runTaskRequest)

    log.info("AWS submission completed:\n{}", runTaskResult)
    val taskArn = runTaskResult.getTasks.asScala.head.getTaskArn

    PendingExecutionHandle(jobDescriptor, AwsJobId(taskArn), None, None)
  }

  private def recover(jobId: BackendJobId): ExecutionHandle = execute()

  private def pollStatus(pendingJob: BackendJobId): BackendRunStatus = {
    val taskArn = pendingJob.taskArn

    val describeTasksRequest = new DescribeTasksRequest()
      .withCluster(clusterName)
      .withTasks(List(taskArn).asJava)

    val describeTasksResult = ecsAsyncClient.describeTasks(describeTasksRequest)

    val tasks = describeTasksResult.getTasks.asScala
    val task = tasks.head
    AwsRunStatus(task)
  }

  private def customPollStatusFailure(oldHandle: BackendPendingExecutionHandle):
  PartialFunction[Exception, ExecutionHandle] = {
    PartialFunction.empty
  }

  private def isTerminal(runStatus: BackendRunStatus): Boolean = {
    runStatus.task.getLastStatus == DesiredStatus.STOPPED.toString
  }

  private def getTerminalMetadata(runStatus: BackendRunStatus): Map[String, Any] = {
    // TODO: Get terminal metadata
    Map.empty
  }

  private def isSuccess(runStatus: BackendRunStatus): Boolean = {
    // TODO: Discriminate failure statuses
    true
  }

  private def handleSuccess(runStatus: BackendRunStatus, handle: BackendPendingExecutionHandle,
                            returnCode: Int): ExecutionHandle = {
    // TODO: Processes success
    log.info("AWS task completed!\n{}", runStatus.task)
    SuccessfulExecutionHandle(Map.empty, returnCode, Map.empty, Seq.empty, None)
  }

  private def handleFailure(runStatus: BackendRunStatus, handle: BackendPendingExecutionHandle): ExecutionHandle = {
    // TODO: Process failure
    FailedNonRetryableExecutionHandle(new Exception(s"Task failed for unknown reason: ${runStatus.task}"), None)
  }

}
