package cromwell.backend.impl.sfs.config

import java.io.PrintWriter
import java.lang.IllegalArgumentException
import java.util.Calendar

import common.validation.Validation._
import cromwell.backend.RuntimeEnvironmentBuilder
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.sfs._
import cromwell.backend.standard.{StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.backend.validation.DockerValidation
import cromwell.core.path.Path
import net.ceedubs.ficus.Ficus._
import wdl.draft2.model._
import wdl.transforms.draft2.wdlom2wom._
import wom.callable.Callable.OptionalInputDefinition
import wom.expression.NoIoFunctionSet
import wom.transforms.WomCommandTaskDefinitionMaker.ops._
import wom.values.{WomEvaluatedCallInputs, WomOptionalValue, WomString, WomValue}

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
    } yield womInputDefinition -> womValue).toMap

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
     *   ${docker} ${script}
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
    if (isDockerRun) {
      Map(
        DockerCwdInput -> WomString(jobPathsWithDocker.callDockerRoot.pathAsString),
        DockerCidInput -> dockerCidInputValue,
        StdoutInput -> WomString(jobPathsWithDocker.toDockerPath(standardPaths.output).pathAsString),
        StderrInput -> WomString(jobPathsWithDocker.toDockerPath(standardPaths.error).pathAsString),
        ScriptInput -> WomString(jobPathsWithDocker.toDockerPath(jobPaths.script).pathAsString),
        JobShellInput -> WomString(jobShell)
      )
    } else {
      Map(
        StdoutInput -> WomString(standardPaths.output.pathAsString),
        StderrInput -> WomString(standardPaths.error.pathAsString),
        ScriptInput -> WomString(jobPaths.script.pathAsString),
        JobShellInput -> WomString(jobShell)
      )
    }
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

  /**
    * Retrieves the DispatchedConfigJob from the stdout using the job-id-regex from the config.
    *
    * @param exitValue The exit value of the dispatch attempt.
    * @param stdout    The stdout from dispatching the job.
    * @param stderr    The stderr from dispatching the job.
    * @return The wrapped job id.
    */
  override def getJob(exitValue: Int, stdout: Path, stderr: Path): StandardAsyncJob = {
    val jobIdRegex = configurationDescriptor.backendConfig.getString(JobIdRegexConfig).r
    val output = stdout.contentAsString.stripLineEnd
    output match {
      case jobIdRegex(jobId) => StandardAsyncJob(jobId)
      case _ =>
        throw new RuntimeException("Could not find job ID from stdout file. " +
          s"Check the stderr file for possible errors: $stderr")
    }
  }

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
    jobScriptArgs(job, "kill", KillTask)
  }

  protected lazy val exitCodeTimeout: Option[Int] = {
    val timeout = configurationDescriptor.backendConfig.as[Option[Int]](ExitCodeTimeoutConfig)
    timeout.foreach{x => if (x < 0) throw new IllegalArgumentException(s"config value '$ExitCodeTimeoutConfig' must be 0 or higher")}
    timeout
  }

  override def pollStatus(handle: StandardAsyncPendingExecutionHandle): SharedFileSystemRunState = {
    handle.previousState match {
      case None =>
        // Is not set yet the status will be set default to running
        SharedFileSystemRunState("Running")
      case Some(s) if (s.status == "Running" || s.status == "WaitingForReturnCode") && jobPaths.returnCode.exists =>
        // If exitcode file does exists status will be set to Done always
        SharedFileSystemRunState("Done")
      case Some(s) if s.status == "Running" =>
        // Exitcode file does not exist at this point, checking is jobs is still alive
        if (exitCodeTimeout.isEmpty) s
        else if (isAlive(handle.pendingJob).fold({ e =>
          log.error(e, s"Running '${checkAliveArgs(handle.pendingJob).argv.mkString(" ")}' did fail")
          true
        }, x => x)) s
        else SharedFileSystemRunState("WaitingForReturnCode")
      case Some(s) if s.status == "WaitingForReturnCode" =>
        // Can only enter this state when the exit code does not exist and the job is not alive anymore
        // `isAlive` is not called anymore from this point

        // If exit-code-timeout is set in the config cromwell will create a fake exitcode file with exitcode 137
        exitCodeTimeout match {
          case Some(timeout) =>
            val currentDate = Calendar.getInstance()
            currentDate.add(Calendar.SECOND, -timeout)
            if (s.date.after(currentDate)) s
            else {
              jobLogger.error(s"Return file not found after $timeout seconds, assuming external kill")
              val writer = new PrintWriter(jobPaths.returnCode.toFile)
              // 137 does mean a external kill -9, this is a assumption but easy workaround for now
              writer.println(9)
              writer.close()
              SharedFileSystemRunState("Failed")
            }
          case _ => s
        }
      case Some(s) if s.status == "Done" => s // Nothing to be done here
      case _ => throw new NotImplementedError("This should not happen, please report this")
    }
  }

  override def isTerminal(runStatus: StandardAsyncRunState): Boolean = {
    runStatus.status == "Done" || runStatus.status == "Failed"
  }

}
