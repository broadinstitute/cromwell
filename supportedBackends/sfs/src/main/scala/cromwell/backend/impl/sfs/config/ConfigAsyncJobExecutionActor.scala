package cromwell.backend.impl.sfs.config

import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.sfs._
import cromwell.backend.standard.{StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.backend.validation.DockerValidation
import cromwell.core.path.Path
import wdl4s.wdl._
import wdl4s.wdl.expression.NoFunctions
import wdl4s.wdl.values.WdlString

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
    val inputsWithFqns = inputs map { case (k, v) => s"$taskName.$k" -> v }
    val command = task.instantiateCommand(task.inputsFromMap(inputsWithFqns), NoFunctions).get
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
      JobNameInput -> WdlString(jobName),
      CwdInput -> WdlString(jobPaths.callRoot.pathAsString)
    )
  }

  /**
    * The inputs that are not specified by the config, that will be passed into a command for either submit or
    * submit-docker.
    */
  private lazy val dockerizedInputs: WorkflowCoercedInputs = {
    if (isDockerRun) {
      Map(
        DockerCwdInput -> WdlString(jobPathsWithDocker.callDockerRoot.pathAsString),
        StdoutInput -> WdlString(jobPathsWithDocker.toDockerPath(jobPaths.stdout).pathAsString),
        StderrInput -> WdlString(jobPathsWithDocker.toDockerPath(jobPaths.stderr).pathAsString),
        ScriptInput -> WdlString(jobPathsWithDocker.toDockerPath(jobPaths.script).pathAsString)
      )
    } else {
      Map(
        StdoutInput -> WdlString(jobPaths.stdout.pathAsString),
        StderrInput -> WdlString(jobPaths.stderr.pathAsString),
        ScriptInput -> WdlString(jobPaths.script.pathAsString)
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
        Option(declarationValidation.key -> WdlString(dockerHash))
      case declarationValidation =>
        declarationValidation.extractWdlValueOption(validatedRuntimeAttributes) map { wdlValue =>
          declarationValidation.key -> wdlValue
        }
    }
    inputOptions.flatten.toMap
  }

  // `runtimeAttributeInputs` has already adjusted for the case of a `JobDescriptor` with `DockerWithHash`.
  override lazy val dockerImageUsed: Option[String] = runtimeAttributeInputs.get(DockerValidation.instance.key).map(_.valueString)
}

/**
  * Submits a job and sends it to the background via "&". Saves the unix PID for status or killing later.
  *
  * @param standardParams Params for running a shared file system job.
  */
class BackgroundConfigAsyncJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends ConfigAsyncJobExecutionActor with BackgroundAsyncJobExecutionActor

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

  /**
    * Generates a command for a job id, using a config task.
    *
    * @param job    The job id.
    * @param suffix The suffix for the scripts.
    * @param task   The config task that defines the command.
    * @return A runnable command.
    */
  private def jobScriptArgs(job: StandardAsyncJob, suffix: String, task: String): SharedFileSystemCommand = {
    val script = jobPaths.script.plusExt(suffix)
    writeTaskScript(script, task, Map(JobIdInput -> WdlString(job.jobId)))
    SharedFileSystemCommand("/bin/bash", script)
  }
}
