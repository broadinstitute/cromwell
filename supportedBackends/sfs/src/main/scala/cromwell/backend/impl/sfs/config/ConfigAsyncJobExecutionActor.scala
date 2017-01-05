package cromwell.backend.impl.sfs.config

import java.nio.file.Path

import better.files._
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.sfs._
import cromwell.backend.standard.{StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.path.PathFactory._
import wdl4s._
import wdl4s.expression.NoFunctions
import wdl4s.values.WdlString

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
    val submitScript = pathPlusSuffix(jobPaths.script, "submit")
    val submitInputs = standardInputs ++ dockerInputs ++ runtimeAttributeInputs
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
  def writeTaskScript(script: File, taskName: String, inputs: WorkflowCoercedInputs): Unit = {
    val task = configInitializationData.wdlNamespace.findTask(taskName).
      getOrElse(throw new RuntimeException(s"Unable to find task $taskName"))
    val inputsWithFqns = inputs map { case (k, v) => s"$taskName.$k" -> v }
    val command = task.instantiateCommand(task.inputsFromMap(inputsWithFqns), NoFunctions).get
    jobLogger.info(s"executing: $command")
    val scriptBody =
      s"""

#!/bin/bash
$command

""".trim + "\n"
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
      CwdInput -> WdlString(jobPaths.callRoot.toString),
      StdoutInput -> WdlString(jobPaths.stdout.toString),
      StderrInput -> WdlString(jobPaths.stderr.toString),
      ScriptInput -> WdlString(jobPaths.script.toString)
    )
  }

  /**
    * Extra arguments if this is a submit-docker command, or Map.empty.
    */
  private lazy val dockerInputs: WorkflowCoercedInputs = {
    if (isDockerRun) {
      Map(
        DockerCwdInput -> WdlString(jobPaths.callDockerRoot.toString)
      )
    } else {
      Map.empty
    }
  }

  /**
    * The arguments generated from the backend config's list of attributes. These will include things like CPU, memory,
    * and other custom arguments like "backend_queue_name", "backend_billing_project", etc.
    */
  private lazy val runtimeAttributeInputs: WorkflowCoercedInputs = {
    val declarationValidations = configInitializationData.declarationValidations
    val inputOptions = declarationValidations map { declarationValidation =>
      declarationValidation.extractWdlValueOption(validatedRuntimeAttributes) map { wdlValue =>
        declarationValidation.key -> wdlValue
      }
    }
    inputOptions.flatten.toMap
  }
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
    val output = File(stdout).contentAsString.stripLineEnd
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
    val script = pathPlusSuffix(jobPaths.script, suffix)
    writeTaskScript(script, task, Map(JobIdInput -> WdlString(job.jobId)))
    SharedFileSystemCommand("/bin/bash", script)
  }
}
