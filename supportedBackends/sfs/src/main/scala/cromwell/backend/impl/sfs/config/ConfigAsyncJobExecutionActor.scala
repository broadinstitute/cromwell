package cromwell.backend.impl.sfs.config

import java.nio.file.FileAlreadyExistsException
import java.time.Instant

import common.validation.Validation._
import cromwell.backend.RuntimeEnvironmentBuilder
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.sfs._
import cromwell.backend.standard.{StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.backend.validation.DockerValidation
import cromwell.core.path.Path
import mouse.all._
import net.ceedubs.ficus.Ficus._
import wdl.draft2.model._
import wdl.transforms.draft2.wdlom2wom._
import wom.callable.Callable.{InputDefinition, OptionalInputDefinition}
import wom.expression.NoIoFunctionSet
import wom.transforms.WomCommandTaskDefinitionMaker.ops._
import wom.values.{WomEvaluatedCallInputs, WomOptionalValue, WomString, WomValue}

import scala.util.{Failure, Success}

/**
  * Base ConfigAsyncJobExecutionActor that reads the config and generates an outer script to submit an inner script
  * containing the call command.
  *
  * There are two types of config executors, 1) those that will require being sent to the unix background, where their
  * PID will be recorded, and 2) those that are dispatched to another system such as Grid Engine, LSF, PBS, etc.
  */
sealed trait ConfigAsyncJobExecutionActor extends SharedFileSystemAsyncJobExecutionActor {

  lazy val configInitializationData: ConfigInitializationData = backendInitializationDataAs[ConfigInitializationData]

  /**
    * Returns the arguments for submitting the job, either with or without docker.
    *
    * At this point, the script might generate a command that isn't dispatched, and it will be up to the
    * BackgroundConfigAsyncJobExecutionActor implementation to send this submission to the background and record the
    * PID. Otherwise, if the submission is dispatched, the other implementation, DispatchedConfigAsyncJobExecutionActor,
    * will grab the job id from the stdout.
    */
  override lazy val processArgs: SharedFileSystemCommand = {
    val submitScript = jobPaths.script.plusExt("submit")
    val submitInputs = standardInputs ++ dockerizedInputs ++ runtimeAttributeInputs
    val submitTaskName = if (isDockerRun) SubmitDockerTask else SubmitTask
    writeTaskScript(submitScript, submitTaskName, submitInputs)
    SharedFileSystemCommand("/bin/bash", submitScript)
  }

  /**
    * Writes a config based task into a file for executing. This might be submit, kill, or even checking if the process
    * is still alive for recover.
    *
    * @param script   Path to write the script.
    * @param taskName The name of the task to retrieve from the precomputed wdl namespace.
    * @param inputs   The customized inputs to this task.
    */
  def writeTaskScript(script: Path, taskName: String, inputs: WorkflowCoercedInputs): Unit = {
    val task = configInitializationData.wdlNamespace.findTask(taskName).
      getOrElse(throw new RuntimeException(s"Unable to find task $taskName"))

    val taskDefinition = task.toWomTaskDefinition.toTry.get

    // Build WOM versions of the provided inputs by correlating with the `InputDefinition`s on the `TaskDefinition`.
    val providedWomInputs: WomEvaluatedCallInputs = (for {
      input <- inputs.toList
      (localName, womValue) = input
      womInputDefinition <- taskDefinition.inputs
      if localName == womInputDefinition.localName.value
    } yield womInputDefinition.asInstanceOf[InputDefinition] -> womValue).toMap

    /*
     * TODO WOM FIXME #2946
     * This forces `none`s onto potentially initialized optionals. This is currently required because otherwise
     * expressions using uninitialized optionals explode on evaluation. e.g.:
     *
     * runtime-attributes = """
     * String? docker
     * String? docker_user
     * """
     * submit = "/bin/bash ${script}"
     * submit-docker = """
     * docker run \
     *   --cidfile ${docker_cid} \
     *   --rm -i \
     *   ${"--user " + docker_user} \
     *   --entrypoint ${job_shell} \
     *   -v ${cwd}:${docker_cwd} \
     *   ${docker} ${docker_script}
     *  """
     *
     * The `docker_user` variable is an uninitialized optional. But the expression `${"--user " + docker_user}` will
     * fail to evaluate if `docker_user` is not provided as an input by the caller of `TaskDefinition#instantiateCommand`.
     * This only appears to be the case when calling `instantiateCommand` on a `TaskDefinition` that came from
     * `WdlTask#womDefinition`, which makes it seem that there is some special initialization of task optionals taking
     * place when create a `WomExecutable` from a `WdlWorkflow`.
     *
     * */
    val optionalsForciblyInitializedToNone = for {
      optional <- taskDefinition.inputs.collect { case oid: OptionalInputDefinition => oid }
      // Only default to `none` if there was no input value explicitly provided.
      if !inputs.contains(optional.localName.value)
    } yield optional -> WomOptionalValue.none(optional.womType.memberType)


    val runtimeEnvironment = RuntimeEnvironmentBuilder(jobDescriptor.runtimeAttributes, jobPaths)(standardParams.minimumRuntimeSettings)
    val allInputs = providedWomInputs ++ optionalsForciblyInitializedToNone
    val womInstantiation = taskDefinition.instantiateCommand(allInputs, NoIoFunctionSet, identity, runtimeEnvironment)

    val command = womInstantiation.toTry.get.commandString
    jobLogger.info(s"executing: $command")
    val scriptBody =
      s"""|#!/bin/bash
          |SCRIPT_COMMAND
          |""".stripMargin.replace("SCRIPT_COMMAND", command)
    script.write(scriptBody)
    ()
  }

  /**
    * The inputs that are not specified by the config, that will be passed into a command for both submit and
    * submit-docker.
    */
  private lazy val standardInputs: WorkflowCoercedInputs = {
    Map(
      JobNameInput -> WomString(jobName),
      CwdInput -> WomString(jobPaths.callRoot.pathAsString)
    )
  }

  private [config] def dockerCidInputValue: WomString = WomString(jobPaths.callExecutionRoot.resolve(jobPaths.dockerCid).pathAsString)

  /**
    * The inputs that are not specified by the config, that will be passed into a command for either submit or
    * submit-docker.
    */
  private lazy val dockerizedInputs: WorkflowCoercedInputs = {
    val dockerPaths = isDockerRun.fold(
      Map(
        DockerCwdInput -> WomString(jobPathsWithDocker.callDockerRoot.pathAsString),
        DockerStdoutInput -> WomString(jobPathsWithDocker.toDockerPath(standardPaths.output).pathAsString),
        DockerStderrInput -> WomString(jobPathsWithDocker.toDockerPath(standardPaths.error).pathAsString),
        DockerScriptInput -> WomString(jobPathsWithDocker.toDockerPath(jobPaths.script).pathAsString),
        DockerCidInput -> dockerCidInputValue,
      ),
      Map.empty
    )

    dockerPaths ++ Map(
      StdoutInput -> WomString(standardPaths.output.pathAsString),
      StderrInput -> WomString(standardPaths.error.pathAsString),
      ScriptInput -> WomString(jobPaths.script.pathAsString),
      JobShellInput -> WomString(jobShell),
    )
  }

  /**
    * The arguments generated from the backend config's list of attributes. These will include things like CPU, memory,
    * and other custom arguments like "backend_queue_name", "backend_billing_project", etc.
    */
  private lazy val runtimeAttributeInputs: WorkflowCoercedInputs = {
    val declarationValidations = configInitializationData.declarationValidations
    val inputOptions = declarationValidations map {
      // Is it always the right thing to pass the Docker hash to a config backend?  What if it can't use hashes?
      case declarationValidation if declarationValidation.key == DockerValidation.instance.key && jobDescriptor.maybeCallCachingEligible.dockerHash.isDefined =>
        val dockerHash = jobDescriptor.maybeCallCachingEligible.dockerHash.get
        Option(declarationValidation.key -> WomString(dockerHash))
      case declarationValidation =>
        declarationValidation.extractWdlValueOption(validatedRuntimeAttributes) map { womValue =>
          declarationValidation.key -> womValue
        }
    }
    inputOptions.flatten.toMap
  }

  // `runtimeAttributeInputs` has already adjusted for the case of a `JobDescriptor` with `DockerWithHash`.
  override lazy val dockerImageUsed: Option[String] = runtimeAttributeInputs.get(DockerValidation.instance.key).map(_.valueString)

  /**
    * Generates a command for a job id, using a config task.
    *
    * @param job    The job id.
    * @param suffix The suffix for the scripts.
    * @param task   The config task that defines the command.
    * @return A runnable command.
    */
  protected def jobScriptArgs(job: StandardAsyncJob, suffix: String, task: String, extraInputs: Map[String, WomValue] = Map.empty): SharedFileSystemCommand = {
    val script = jobPaths.script.plusExt(suffix)
    writeTaskScript(script, task, Map(JobIdInput -> WomString(job.jobId)) ++ extraInputs)
    SharedFileSystemCommand("/bin/bash", script)
  }
}

/**
  * Submits a job and sends it to the background via "&". Saves the unix PID for status or killing later.
  *
  * @param standardParams Params for running a shared file system job.
  */
class BackgroundConfigAsyncJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends ConfigAsyncJobExecutionActor with BackgroundAsyncJobExecutionActor {

  override def killArgs(job: StandardAsyncJob): SharedFileSystemCommand = {
    if (isDockerRun) jobScriptArgs(job, "kill", KillDockerTask, Map(DockerCidInput -> dockerCidInputValue))
    else super[BackgroundAsyncJobExecutionActor].killArgs(job)
  }
}

/**
  * Submits a job and returns relatively quickly. The job-id-regex is then used to read the job id for status or killing
  * later.
  *
  * @param standardParams Params for running a shared file system job.
  */
class DispatchedConfigAsyncJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends ConfigAsyncJobExecutionActor {

  lazy val jobIdRegexString = configurationDescriptor.backendConfig.getString(JobIdRegexConfig)

  /**
    * Retrieves the DispatchedConfigJob from the stdout using the job-id-regex from the config.
    *
    * @param exitValue The exit value of the dispatch attempt.
    * @param stdout    The stdout from dispatching the job.
    * @param stderr    The stderr from dispatching the job.
    * @return The wrapped job id.
    */
  override def getJob(exitValue: Int, stdout: Path, stderr: Path): StandardAsyncJob =
    DispatchedConfigAsyncJobExecutionActor.getJob(stdout.contentAsString, stderr, jobIdRegexString)

  /**
    * Checks if the job is alive using the command from the config.
    *
    * @param job The job to check.
    * @return A command that checks if the job is alive.
    */
  override def checkAliveArgs(job: StandardAsyncJob): SharedFileSystemCommand = {
    jobScriptArgs(job, "check", CheckAliveTask)
  }

  /**
    * Kills the job using the kill command from the config.
    *
    * @param job The job id to kill.
    * @return A command that may be used to kill the job.
    */
  override def killArgs(job: StandardAsyncJob): SharedFileSystemCommand = {
    if (isDockerRun) jobScriptArgs(job, "kill", KillDockerTask, Map(DockerCidInput -> dockerCidInputValue))
    else jobScriptArgs(job, "kill", KillTask)
  }

  protected lazy val exitCodeTimeout: Option[Long] = {
    val timeout = configurationDescriptor.backendConfig.as[Option[Long]](ExitCodeTimeoutConfig)
    timeout match {
      case Some(x) =>
        jobLogger.info("Cromwell will watch for an rc file *and* double-check every {} seconds to make sure this job is still alive", x)
        if (x < 0) throw new IllegalArgumentException(s"config value '$ExitCodeTimeoutConfig' must be 0 or higher")
      case None =>
        jobLogger.info("Cromwell will watch for an rc file but will *not* double-check whether this job is actually alive (unless Cromwell restarts)")
    }
    timeout
  }

  override def pollStatus(handle: StandardAsyncPendingExecutionHandle): SharedFileSystemRunState = {

    lazy val nextTimeout = exitCodeTimeout.map(t => Instant.now.plusSeconds(t))

    handle.previousState match {
      case None =>
        // Is not set yet the status will be set default to running
        SharedFileSystemJobRunning(nextTimeout)

      case Some(_) if jobPaths.returnCode.exists =>
        // If exitcode file does exists status will be set to Done always
        SharedFileSystemJobDone

      case Some(running: SharedFileSystemJobRunning) =>
        if (running.stale) {
          // Since the exit code timeout has expired for a running job, check whether the job is still alive.

          isAlive(handle.pendingJob) match {
            case Success(true) =>
              // The job is still running. Don't check again for another timeout.
              SharedFileSystemJobRunning(nextTimeout)
            case Success(false) =>
              // The job has stopped but we don't have an RC yet. We'll wait one more 'timeout' for the RC to arrive:
              SharedFileSystemJobWaitingForReturnCode(nextTimeout)
            case Failure(e) =>
              log.error(e, s"Failed to check status for ${handle.jobDescriptor.key.tag} using command: ${checkAliveArgs(handle.pendingJob)}")
              SharedFileSystemJobRunning(nextTimeout)
          }
        } else {
          // Not stale yet so keep on running!
          running
        }

      case Some(waiting: SharedFileSystemJobWaitingForReturnCode) =>
        if (waiting.giveUpWaiting) {
          // Can only enter this state when the exit code does not exist and the job is not alive anymore
          // `isAlive` is not called anymore from this point

          // If exit-code-timeout is set in the config cromwell will create a fake exitcode file
          val backupError = "??? (!! Programmer Error: It should be impossible to give up on 'waiting' without having set a maximum wait timeout. Please report this as a bug in the Cromwell Github repository !!)"
          jobLogger.error(s"Return file not found after ${exitCodeTimeout.getOrElse(backupError)} seconds, assuming external kill")

          val returnCodeTemp = jobPaths.returnCode.plusExt("kill")

          // If SIGTERM, SIGKILL or SIGINT codes are used, cromwell will assume the job has been aborted by cromwell.
          // And it will therefore NOT retry the job. Which makes perfect sense. Best not to change that in the
          // StandardAsyncExecutionActor code, but give a user-defined return code here.
          // http://tldp.org/LDP/abs/html/exitcodes.html gives information on exit codes and suggests restricting
          // user-defined exit codes to the range 64-113.
          // Since it is arbitrary which code is chosen from that range, and it has to relate with the unpleasant
          // business of 'killing'. 79 was chosen. The year that video killed the radio star: https://youtu.be/W8r-tXRLazs
          returnCodeTemp.appendLine("79")

          try {
            returnCodeTemp.moveTo(jobPaths.returnCode)
            SharedFileSystemJobFailed
          } catch {
            case _: FileAlreadyExistsException =>
              // If the process has already completed, there will be an existing RC:
              // Essentially, this covers the rare race condition whereby the task completes between starting to write the
              // fake RC, and trying to copy it:

              log.error(s"An RC file appeared at ${jobPaths.returnCode} whilst trying to copy a fake exitcode file from ${returnCodeTemp}. Not to worry: the real file should now be picked up on the next poll.")

              // Delete the fake file since it's not needed:
              returnCodeTemp.delete(true)

              // We shouldn't be here anymore...
              // Return a "no-op" for now and allow the now-existing rc file to be picked up as normal on the next
              // iteration through the polling logic:
              waiting
          }
        } else {
          // Not giving up yet so keep on waiting... for now...
          waiting
        }

      case Some(otherStatus) =>
        jobLogger.warn(s"Unexpected status poll received when already in status $otherStatus")
        otherStatus
    }
  }

  override def isTerminal(runStatus: SharedFileSystemRunState): Boolean = runStatus.terminal
}


object DispatchedConfigAsyncJobExecutionActor {
  def getJob(stdoutContent: String, stderr: Path, jobIdRegexString: String): StandardAsyncJob = {
    val jobIdRegex = jobIdRegexString.r
    val output = stdoutContent.stripLineEnd
    jobIdRegex findFirstIn output match {
      case Some(jobIdRegex(jobId)) => StandardAsyncJob(jobId)
      case _ =>
        throw new RuntimeException("Could not find job ID from stdout file." +
          s"Check the stderr file for possible errors: $stderr")
    }
  }
}

