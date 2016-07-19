package cromwell.backend.impl.config

import java.nio.file.Path

import better.files._
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.sfs._

/**
  * Builds a backend by reading the job control from the config.
  *
  * This class is considered experimental.
  *
  * In the future we may:
  *   - Parse submit and kill command as a wdl expression
  *   - Allow a syntax to specify runtime attribute validations based on a set of standard keys
  *   - Allow a syntax to specify runtime attribute validations based on a regex string
  *
  * @param configurationDescriptor The config information.
  */
class ConfigBackendLifecycleActorFactory(override val configurationDescriptor: BackendConfigurationDescriptor)
  extends SharedFileSystemBackendLifecycleActorFactory {

  {
    val required = Seq("submit", "kill", "job-id-regex")
    val hasAll = required forall configurationDescriptor.backendConfig.hasPath
    if (!hasAll) throw new RuntimeException(s"Config must contain strings for keys '${required.mkString("', '")}'")
  }

  override def asyncJobExecutionActorClass = classOf[ConfigBackendAsyncJobExecutionActor]
}

class ConfigBackendAsyncJobExecutionActor(override val params: SharedFileSystemAsyncJobExecutionActorParams)
  extends SharedFileSystemAsyncJobExecutionActor[ConfigBackendJob] {

  def configString(key: String) = configurationDescriptor.backendConfig.getString(key)

  def configArgs(key: String) = configString(key).split("\\s+")

  // Primitive check for a wdl type key
  def isKey(arg: String, key: String) = arg == "${%s}".format(key)

  override def processArgs = {
    val argv = configArgs("submit") map {
      case arg if isKey(arg, "job_name") => jobName
      case arg if isKey(arg, "cwd") => jobPaths.callRoot.toAbsolutePath
      case arg if isKey(arg, "out") => jobPaths.stdout.toAbsolutePath
      case arg if isKey(arg, "err") => jobPaths.stderr.toAbsolutePath
      case arg if isKey(arg, "script") => jobPaths.script.toAbsolutePath
      case arg => arg
    }
    new SharedFileSystemCommand(argv)
  }

  override def getJob(exitValue: Int, stdout: Path, stderr: Path) = {
    val jobIdRegex = configString("job-id-regex").r
    val output = stdout.contentAsString.stripLineEnd
    output match {
      case jobIdRegex(jobId) => new ConfigBackendJob(jobId)
      case _ =>
        throw new RuntimeException("Could not find job ID from stdout file. " +
          s"Check the stderr file for possible errors: ${stderr.toAbsolutePath}")
    }
  }

  override def killArgs(job: ConfigBackendJob) = {
    val argv = configArgs("kill") map {
      case arg if isKey(arg, "job_id") => job.jobId
      case arg => arg
    }

    new SharedFileSystemCommand(argv)
  }
}

class ConfigBackendJob(override val jobId: String) extends SharedFileSystemJob
