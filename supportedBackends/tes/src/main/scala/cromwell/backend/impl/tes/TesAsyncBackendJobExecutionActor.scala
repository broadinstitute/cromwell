package cromwell.backend.impl.tes

import akka.http.scaladsl.Http
import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport._
import akka.http.scaladsl.marshalling.Marshal
import akka.http.scaladsl.model.HttpHeader.ParsingResult.Ok
import akka.http.scaladsl.model._
import akka.http.scaladsl.unmarshalling.{Unmarshal, Unmarshaller}
import akka.stream.ActorMaterializer
import akka.util.ByteString
import cats.implicits._
import common.collections.EnhancedCollections._
import common.exception.AggregatedMessageException
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cromwell.backend.BackendJobLifecycleActor
import cromwell.backend.async.{
  AbortedExecutionHandle,
  ExecutionHandle,
  FailedNonRetryableExecutionHandle,
  PendingExecutionHandle
}
import cromwell.backend.impl.tes.TesAsyncBackendJobExecutionActor.{
  determineWSMSasEndpointFromInputs,
  generateLocalizedSasScriptPreamble
}
import cromwell.backend.impl.tes.TesResponseJsonFormatter._
import cromwell.backend.standard.{
  ScriptPreambleData,
  StandardAsyncExecutionActor,
  StandardAsyncExecutionActorParams,
  StandardAsyncJob
}
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.retry.Retry._
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.filesystems.blob.{BlobPath, WSMBlobSasTokenGenerator}
import cromwell.filesystems.drs.{DrsPath, DrsResolver}
import net.ceedubs.ficus.Ficus._
import wom.values.WomFile

import java.io.FileNotFoundException
import java.nio.file.FileAlreadyExistsException
import java.time.Duration
import java.time.temporal.ChronoUnit
import scala.concurrent.Future
import scala.util.{Failure, Success, Try}
sealed trait TesRunStatus {
  def isTerminal: Boolean
  def sysLogs: Seq[String] = Seq.empty[String]
}

case object Running extends TesRunStatus {
  def isTerminal = false
}

case object Complete extends TesRunStatus {
  def isTerminal = true
}

case class Error(override val sysLogs: Seq[String] = Seq.empty[String]) extends TesRunStatus {
  def isTerminal = true
  override def toString = "SYSTEM_ERROR"
}

case class Failed(override val sysLogs: Seq[String] = Seq.empty[String]) extends TesRunStatus {
  def isTerminal = true
  override def toString = "EXECUTOR_ERROR"
}

case object Cancelled extends TesRunStatus {
  def isTerminal = true
}

object TesAsyncBackendJobExecutionActor {
  val JobIdKey = "tes_job_id"

  def generateLocalizedSasScriptPreamble(environmentVariableName: String, getSasWsmEndpoint: String): String =
    // BEARER_TOKEN: https://learn.microsoft.com/en-us/azure/active-directory/managed-identities-azure-resources/how-to-use-vm-token#get-a-token-using-http
    // NB: Scala string interpolation and bash variable substitution use similar syntax. $$ is an escaped $, much like \\ is an escaped \.
    s"""
       |### BEGIN ACQUIRE LOCAL SAS TOKEN ###
       |# Function to check if a command exists on this machine
       |command_exists() {
       |  command -v "$$1" > /dev/null 2>&1
       |}
       |
       |# Require the user image to pre-install `jq` and `curl` for us. Empirically, we found
       |# that `apt install` at scale can run into repo outages and flakiness.
       |if ! command_exists curl; then
       |  echo "Error: user image must include `curl` for just-in-time SAS token generation, but it is not installed."
       |  exit 1
       |fi
       |
       |if ! command_exists jq; then
       |  echo "Error: user image must include `jq` for just-in-time SAS token generation, but it is not installed."
       |  exit 1
       |fi
       |curl --version
       |jq --version
       |# Acquire bearer token, relying on the User Assigned Managed Identity of this VM.
       |echo Acquiring Bearer Token using User Assigned Managed Identity...
       |BEARER_TOKEN=$$(curl 'http://169.254.169.254/metadata/identity/oauth2/token?api-version=2018-02-01&resource=https%3A%2F%2Fmanagement.azure.com%2F' -H Metadata:true -s | jq .access_token)
       |# Remove the leading and trailing quotes
       |BEARER_TOKEN="$${BEARER_TOKEN#\\"}"
       |BEARER_TOKEN="$${BEARER_TOKEN%\\"}"
       |
       |# Use the precomputed endpoint from cromwell + WSM to acquire a sas token
       |echo Requesting sas token from WSM...
       |sas_response_json=$$(curl -s \\
       |                    --retry 3 \\
       |                    --retry-delay 2 \\
       |                    -X POST "$getSasWsmEndpoint" \\
       |                    -H "Content-Type: application/json" \\
       |                    -H "accept: */*" \\
       |                    -H "Authorization: Bearer $${BEARER_TOKEN}" \\
       |                    -H "Content-Length: 0" \\
       |                    -d "")
       |
       |# Store token as environment variable
       |export $environmentVariableName=$$(echo "$${sas_response_json}" | jq -r '.token')
       |
       |# Echo the first characters for logging/debugging purposes. "null" indicates something went wrong.
       |echo "Saving sas token: $${$environmentVariableName:0:4}**** to environment variable $environmentVariableName..."
       |### END ACQUIRE LOCAL SAS TOKEN ###
       |""".stripMargin

  private def maybeConvertToBlob(pathToTest: Try[Path]): Try[BlobPath] =
    pathToTest.collect { case blob: BlobPath => blob }

  /**
   * Computes an endpoint that can be used to retrieve a sas token for a particular blob storage container.
   * This assumes that some of the task inputs are blob files, all blob files are in the same container, and we can get a sas
   * token for this container from WSM.
   * The task VM will use the user assigned managed identity that it is running as in order to authenticate.
   * @param taskInputs The inputs to this particular TesTask. If any are blob files, the first  will be used to
   *                   determine the storage container to retrieve the sas token for.
   * @param pathGetter A function to convert string filepath into a cromwell Path object.
   * @param blobConverter A function to convert a Path into a Blob path, if possible. Provided for testing purposes.
   * @return A URL endpoint that, when called with proper authentication, will return a sas token.
   *         Returns 'None' if one should not be used for this task.
   */
  def determineWSMSasEndpointFromInputs(taskInputs: List[Input],
                                        pathGetter: String => Try[Path],
                                        logger: JobLogger,
                                        blobConverter: Try[Path] => Try[BlobPath] = maybeConvertToBlob
  ): Try[String] = {
    // Collect all of the inputs that are valid blob paths
    val blobFiles = taskInputs
      .collect { case Input(_, _, Some(url), _, _, _) => blobConverter(pathGetter(url)) }
      .collect { case Success(blob) => blob }

    // Log if not all input files live in the same container.
    if (blobFiles.map(_.container).distinct.size > 1) {
      logger.info(
        s"While parsing blob inputs, found more than one container. Generating SAS token based on first file in the list."
      )
    }

    // We use the first blob file in the list to determine the correct blob container.
    blobFiles.headOption
      .map { blobPath =>
        blobPath.getFilesystemManager.blobTokenGenerator match {
          case wsmGenerator: WSMBlobSasTokenGenerator =>
            wsmGenerator.getWSMSasFetchEndpoint(blobPath, Some(Duration.of(24, ChronoUnit.HOURS)))
          case _ =>
            Failure(new UnsupportedOperationException("Blob file does not have an associated WSMBlobSasTokenGenerator"))
        }
      }
      .getOrElse(
        Failure(
          new NoSuchElementException(
            "Could not infer blob storage container from task inputs: No valid blob files provided."
          )
        )
      )
  }
}

class TesAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
    extends BackendJobLifecycleActor
    with StandardAsyncExecutionActor
    with TesJobCachingActorHelper {
  implicit val actorSystem = context.system
  implicit val materializer = ActorMaterializer()

  override type StandardAsyncRunInfo = Any

  override type StandardAsyncRunState = TesRunStatus

  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that
  override lazy val pollBackOff: SimpleExponentialBackoff = tesConfiguration.pollBackoff
  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = tesConfiguration.executeOrRecoverBackoff

  private lazy val realDockerImageUsed: String =
    jobDescriptor.maybeCallCachingEligible.dockerHash.getOrElse(runtimeAttributes.dockerImage)
  override lazy val dockerImageUsed: Option[String] = Option(realDockerImageUsed)

  private val tesEndpoint = workflowDescriptor.workflowOptions.getOrElse("endpoint", tesConfiguration.endpointURL)

  override lazy val jobTag: String = jobDescriptor.key.tag

  private val outputMode = validate {
    OutputMode.withName(
      configurationDescriptor.backendConfig
        .getAs[String]("output-mode")
        .getOrElse("granular")
        .toUpperCase
    )
  }

  /**
   * This script preamble is bash code that is executed at the start of a task inside the user's container.
   * It is executed directly before the user's instantiated command is, which gives cromwell a chance to adjust the
   * container environment before the actual task runs. See commandScriptContents in StandardAsyncExecutionActor for more context.
   *
   * For TES tasks, we sometimes want to acquire and save an azure sas token to an environment variable.
   * If the user provides a value for runtimeAttributes.localizedSasEnvVar, we will add the relevant bash code to the preamble
   * that acquires/exports the sas token to an environment variable. Once there, it will be visible to the user's task code.
   *
   * If runtimeAttributes.localizedSasEnvVar is provided in the WDL (and determineWSMSasEndpointFromInputs is successful),
   * we will export the sas token to an environment variable named to be the value of runtimeAttributes.localizedSasEnvVar.
   * Otherwise, we won't alter the preamble.
   *
   * See determineWSMSasEndpointFromInputs to see how we use taskInputs to infer *which* container to get a sas token for.
   *
   * @return Bash code to run at the start of a task.
   */
  override def scriptPreamble: ErrorOr[ScriptPreambleData] =
    runtimeAttributes.localizedSasEnvVar match {
      case Some(environmentVariableName) =>
        { // Case: user wants a sas token. Return the computed preamble or die trying.
          val workflowName = workflowDescriptor.callable.name
          val callInputFiles = jobDescriptor.fullyQualifiedInputs.safeMapValues {
            _.collectAsSeq { case w: WomFile => w }
          }
          val taskInputs: List[Input] = TesTask.buildTaskInputs(callInputFiles, workflowName, mapCommandLineWomFile)
          val computedEndpoint = determineWSMSasEndpointFromInputs(taskInputs, getPath, jobLogger)
          computedEndpoint.map(endpoint =>
            ScriptPreambleData(generateLocalizedSasScriptPreamble(environmentVariableName, endpoint),
                               executeInSubshell = false
            )
          )
        }.toErrorOr
      case _ =>
        ScriptPreambleData("",
                           executeInSubshell = false
        ).valid // Case: user doesn't want a sas token. Empty preamble is the correct preamble.
    }
  override def mapCommandLineWomFile(womFile: WomFile): WomFile =
    womFile.mapFile(value =>
      (getPath(value), asAdHocFile(womFile)) match {
        case (Success(path: Path), Some(adHocFile)) =>
          // Ad hoc files will be placed directly at the root ("/cromwell_root/ad_hoc_file.txt") unlike other input files
          // for which the full path is being propagated ("/cromwell_root/path/to/input_file.txt")
          tesJobPaths.containerExec(commandDirectory, adHocFile.alternativeName.getOrElse(path.name))
        case _ => mapCommandLineJobInputWomFile(womFile).value
      }
    )

  override def mapCommandLineJobInputWomFile(womFile: WomFile): WomFile =
    womFile.mapFile(value =>
      getPath(value) match {
        case Success(drsPath: DrsPath) =>
          val filepath = DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
          tesJobPaths.containerExec(commandDirectory, filepath)
        case Success(path: Path) if path.startsWith(tesJobPaths.workflowPaths.DockerRoot) =>
          path.pathAsString
        case Success(path: Path) if path.equals(tesJobPaths.callExecutionRoot) =>
          commandDirectory.pathAsString
        case Success(path: Path) if path.startsWith(tesJobPaths.callExecutionRoot) =>
          tesJobPaths.containerExec(commandDirectory, path.name)
        case Success(path: Path) if path.startsWith(tesJobPaths.callRoot) =>
          tesJobPaths.callDockerRoot.resolve(path.name).pathAsString
        case Success(path: BlobPath) if path.startsWith(tesJobPaths.workflowPaths.workflowRoot) =>
          // Blob paths can get really long, which causes problems for some tools. If this input file
          // lives in the workflow execution directory, strip off that prefix from the path we're
          // generating inside `inputs/` to keep the total path length under control.
          // In Terra on Azure, this saves us 200+ characters.
          tesJobPaths.callInputsDockerRoot
            .resolve(
              path.pathStringWithoutPrefix(tesJobPaths.workflowPaths.workflowRoot)
            )
            .pathAsString
        case Success(path: BlobPath) if path.startsWith(tesJobPaths.workflowPaths.executionRoot) =>
          // See comment above... if this file is in the execution root, strip that off.
          // In Terra on Azure, this saves us 160+ characters.
          tesJobPaths.callInputsDockerRoot
            .resolve(
              path.pathStringWithoutPrefix(tesJobPaths.workflowPaths.executionRoot)
            )
            .pathAsString
        case Success(path: Path) =>
          tesJobPaths.callInputsDockerRoot.resolve(path.pathWithoutScheme.stripPrefix("/")).pathAsString
        case _ =>
          value
      }
    )

  override lazy val commandDirectory: Path =
    runtimeAttributes.dockerWorkingDir match {
      case Some(path) => DefaultPathBuilder.get(path)
      case None => tesJobPaths.callExecutionDockerRoot
    }

  def createTaskMessage(): ErrorOr[Task] = {
    val tesTask = (commandScriptContents, outputMode).mapN { case (contents, mode) =>
      TesTask(
        jobDescriptor,
        configurationDescriptor,
        jobLogger,
        tesJobPaths,
        runtimeAttributes,
        commandDirectory,
        contents,
        instantiatedCommand,
        realDockerImageUsed,
        mapCommandLineWomFile,
        jobShell,
        mode
      )
    }

    tesTask.map(TesTask.makeTask)
  }

  def writeScriptFile(): Future[Unit] =
    commandScriptContents.fold(
      errors => Future.failed(new RuntimeException(errors.toList.mkString(", "))),
      asyncIo.writeAsync(jobPaths.script, _, Seq.empty)
    )

  override def executeAsync(): Future[ExecutionHandle] = {
    // create call exec dir
    tesJobPaths.callExecutionRoot.createPermissionedDirectories()
    val taskMessageFuture = createTaskMessage().fold(
      errors => Future.failed(new RuntimeException(errors.toList.mkString(", "))),
      Future.successful
    )

    for {
      _ <- writeScriptFile()
      taskMessage <- taskMessageFuture
      entity <- Marshal(taskMessage).to[RequestEntity]
      ctr <- makeRequest[CreateTaskResponse](HttpRequest(method = HttpMethods.POST, uri = tesEndpoint, entity = entity))
    } yield PendingExecutionHandle(jobDescriptor, StandardAsyncJob(ctr.id), None, previousState = None)
  }

  override def reconnectAsync(jobId: StandardAsyncJob) = {
    val handle = PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState](jobDescriptor,
                                                                                                       jobId,
                                                                                                       None,
                                                                                                       previousState =
                                                                                                         None
    )
    Future.successful(handle)
  }

  override def recoverAsync(jobId: StandardAsyncJob) = reconnectAsync(jobId)

  override def reconnectToAbortAsync(jobId: StandardAsyncJob) = {
    tryAbort(jobId)
    reconnectAsync(jobId)
  }

  override def tryAbort(job: StandardAsyncJob): Unit = {

    val returnCodeTmp = jobPaths.returnCode.plusExt("kill")
    returnCodeTmp.write(s"$SIGTERM\n")
    try
      returnCodeTmp.moveTo(jobPaths.returnCode)
    catch {
      case _: FileAlreadyExistsException =>
        // If the process has already completed, there will be an existing rc file.
        returnCodeTmp.delete(true)
    }

    makeRequest[CancelTaskResponse](
      HttpRequest(method = HttpMethods.POST, uri = s"$tesEndpoint/${job.jobId}:cancel")
    ) onComplete {
      case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, job.jobId)
      case Failure(ex) => jobLogger.warn("{} Failed to abort {}", tag, job.jobId, ex)
    }

    ()
  }

  override def requestsAbortAndDiesImmediately: Boolean = false

  override def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[TesRunStatus] =
    for {
      status <- queryStatusAsync(handle)
      errorLog <- status match {
        case Error(_) | Failed(_) => getErrorLogs(handle)
        case _ => Future.successful(Seq.empty[String])
      }
      statusWithLog = status match {
        case Error(_) => Error(errorLog)
        case Failed(_) => Failed(errorLog)
        case _ => status
      }
    } yield statusWithLog

  private def queryStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[TesRunStatus] =
    makeRequest[MinimalTaskView](HttpRequest(uri = s"$tesEndpoint/${handle.pendingJob.jobId}?view=MINIMAL")) map {
      response =>
        val state = response.state
        state match {
          case s if s.contains("COMPLETE") =>
            jobLogger.info(s"Job ${handle.pendingJob.jobId} is complete")
            Complete

          case s if s.contains("CANCELED") =>
            jobLogger.info(s"Job ${handle.pendingJob.jobId} was canceled")
            Cancelled

          case s if s.contains("EXECUTOR_ERROR") =>
            jobLogger.info(s"TES reported a failure for Job ${handle.pendingJob.jobId}: '$s'")
            Failed()

          case s if s.contains("SYSTEM_ERROR") =>
            jobLogger.info(s"TES reported an error for Job ${handle.pendingJob.jobId}: '$s'")
            Error()

          case _ => Running
        }
    }

  private def getErrorLogs(handle: StandardAsyncPendingExecutionHandle): Future[Seq[String]] =
    makeRequest[Task](HttpRequest(uri = s"$tesEndpoint/${handle.pendingJob.jobId}?view=FULL")) map { response =>
      response.logs.flatMap(_.lastOption).flatMap(_.system_logs).getOrElse(Seq.empty[String])
    }

  override def customPollStatusFailure: PartialFunction[(ExecutionHandle, Exception), ExecutionHandle] = {
    case (oldHandle: StandardAsyncPendingExecutionHandle @unchecked, e: Exception) =>
      jobLogger.error(s"$tag TES Job ${oldHandle.pendingJob.jobId} has not been found, failing call")
      FailedNonRetryableExecutionHandle(e, kvPairsToSave = None)
  }

  private def handleExecutionError(status: TesRunStatus, returnCode: Option[Int]): Future[ExecutionHandle] = {
    val msg = s"Task ${jobDescriptor.key.tag} failed with ${status.toString}"
    jobLogger.info(s"${msg}. Error messages: ${status.sysLogs}")
    val exception = new AggregatedMessageException(msg, status.sysLogs)
    Future.successful(FailedNonRetryableExecutionHandle(exception, returnCode, None))
  }

  override def handleExecutionFailure(status: StandardAsyncRunState, returnCode: Option[Int]) =
    status match {
      case Cancelled => Future.successful(AbortedExecutionHandle)
      case Error(_) | Failed(_) => handleExecutionError(status, returnCode)
      case _ => super.handleExecutionFailure(status, returnCode)
    }

  override def isTerminal(runStatus: TesRunStatus): Boolean =
    runStatus.isTerminal

  override def isDone(runStatus: TesRunStatus): Boolean =
    runStatus match {
      case Complete => true
      case _ => false
    }

  override def mapOutputWomFile(womFile: WomFile): WomFile =
    womFile mapFile { path =>
      val absPath = getPath(path) match {
        case Success(absoluteOutputPath) if absoluteOutputPath.isAbsolute => absoluteOutputPath
        case _ => tesJobPaths.callExecutionRoot.resolve(path)
      }

      if (!absPath.exists) {
        throw new FileNotFoundException(s"Could not process output, file not found: ${absPath.pathAsString}")
      } else absPath.pathAsString
    }

  // Headers that should be included with all requests to the TES server
  private def requestHeaders: List[HttpHeader] =
    tesConfiguration.token.flatMap { t =>
      HttpHeader.parse("Authorization", t) match {
        case Ok(header, _) => Some(header)
        case _ => None
      }
    }.toList

  private def makeRequest[A](request: HttpRequest)(implicit um: Unmarshaller[ResponseEntity, A]): Future[A] =
    for {
      response <- withRetry(() => Http().singleRequest(request.withHeaders(requestHeaders)))
      data <-
        if (response.status.isFailure()) {
          response.entity.dataBytes.runFold(ByteString(""))(_ ++ _).map(_.utf8String) flatMap { errorBody =>
            Future.failed(
              new RuntimeException(s"Failed TES request: Code ${response.status.intValue()}, Body = $errorBody")
            )
          }
        } else {
          Unmarshal(response.entity).to[A]
        }
    } yield data
}
