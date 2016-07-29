package cromwell.backend.impl.sge

import java.nio.file.Path

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.sfs._
import cromwell.backend.validation.MemoryValidation
import wdl4s.parser.MemoryUnit

class SgeBackendLifecycleActorFactory(override val configurationDescriptor: BackendConfigurationDescriptor)
  extends SharedFileSystemBackendLifecycleActorFactory {

  override def asyncJobExecutionActorClass = classOf[SgeAsyncJobExecutionActor]

  override def runtimeAttributesBuilder =
    super.runtimeAttributesBuilder.withValidation(MemoryValidation.optional)
}

class SgeAsyncJobExecutionActor(override val params: SharedFileSystemAsyncJobExecutionActorParams)
  extends SharedFileSystemAsyncJobExecutionActor[SgeJob] {

  override def processArgs: SharedFileSystemCommand = {
    import lenthall.config.ScalaConfig._

    val config = configurationDescriptor.backendConfig
    val queueParam = config.getStringOption("queue").map(Seq("-q", _)).getOrElse(Seq.empty)
    val projectParam = config.getStringOption("project").map(Seq("-P", _)).getOrElse(Seq.empty)
    val memoryParamValue = for {
      memoryName <- config.getStringOption("memoryParam")
      memoryValue <- MemoryValidation.optional.extract(validatedRuntimeAttributes)
    } yield s"$memoryName=${memoryValue.to(MemoryUnit.GB)}G"
    val memoryParam = memoryParamValue.map(Seq("-l", _)).getOrElse(Seq.empty)

    val argv = Seq(
      "qsub",
      "-terse",
      "-N", jobName,
      "-V",
      "-b", "n",
      "-wd", jobPaths.callRoot.toAbsolutePath,
      "-o", jobPaths.stdout.toAbsolutePath,
      "-e", jobPaths.stderr.toAbsolutePath) ++
      queueParam ++
      projectParam ++
      memoryParam ++
      Seq(jobPaths.script.toAbsolutePath)

    new SharedFileSystemCommand(argv)
  }

  def getJob(exitValue: Int, stdout: Path, stderr: Path): SgeJob = {
    import better.files._
    val jobId = stdout.contentAsString.stripLineEnd.toInt
    new SgeJob(jobId.toString)
  }

  override def killArgs(job: SgeJob): SharedFileSystemCommand = {
    SharedFileSystemCommand("qdel", job.jobId)
  }
}

class SgeJob(override val jobId: String) extends SharedFileSystemJob
