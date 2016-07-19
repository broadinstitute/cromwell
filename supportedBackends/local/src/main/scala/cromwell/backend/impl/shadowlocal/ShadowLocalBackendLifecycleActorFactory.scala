package cromwell.backend.impl.shadowlocal

import java.nio.file.Path

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.sfs._
import cromwell.backend.validation.DockerValidation
import cromwell.core.PathFactory._
import cromwell.core.PathWriter

import scala.sys.process.Process

/**
  * An experimental / beta / first pass of moving Local Backend to the Shared File System (SFS).
  *
  * The goal is to have the functionality now in Local Backend as just a thin facade to the SFS. This could either be as
  * the current SGE backend, or hopefully able to even use the Config-based backend- assuming the code is common enough.
  * If we can't completely remove the concrete implementation, the secondary goal is to keep the Local Backend sharing
  * as much code with the SFS.
  *
  * Currently, this version of the Shadow Local Backend runs the process just like the original Local Backend, including
  * attaching threads and draining the stdout/stderr.
  *
  * http://stackoverflow.com/questions/16983372/why-does-process-hang-if-the-parent-does-not-consume-stdout-stderr-in-java
  *
  * After the process is launched, the UNIX PID is recorded using reflection, and the process object otherwise forgotten
  * about by this Shadow Local Backend.
  *
  * When the process needs to be killed, instead of calling `Process.destroy()`, the unix command kill is invoked
  * instead.
  *
  * It should be straight-forward enough to have a concept of an "attached" backend, even down in the SFS. We may
  * need another yet-another wrapper script though.
  *
  * http://stackoverflow.com/questions/931536/how-do-i-launch-a-completely-independent-process-from-a-java-program
  *
  * This would allow any backend that runs "attached", perhaps even an ssh based backend.
  *
  * For now, the SFS doesn't know about the concept of "attached" shared file system backends, so we must override
  * certain methods to keep the base SFS from trying to close the stdout/stderr after process launch. Again, in the
  * future, hopefully the "attached" backends wouldn't attach themselves. A wrapper script could behave like a "local"
  * implementation of "qsub -terse": Redirect a sub-process' stdout / stderr to the specified files, and just return the
  * sub-process PID in a string in the stdout. We wouldn't then need reflection, as we wouldn't care about the wrapper
  * PID.
  *
  * @param configurationDescriptor The backend and global config.
  */
class ShadowLocalBackendLifecycleActorFactory(override val configurationDescriptor: BackendConfigurationDescriptor)
  extends SharedFileSystemBackendLifecycleActorFactory {
  override def asyncJobExecutionActorClass = classOf[ShadowLocalAsyncJobExecutionActor]

  override def supportsDocker = true
}

case class ShadowLocalJob(jobId: String) extends SharedFileSystemJob

class ShadowLocalAsyncJobExecutionActor(override val params: SharedFileSystemAsyncJobExecutionActorParams)
  extends SharedFileSystemAsyncJobExecutionActor[ShadowLocalJob] {

  override def processArgs = {
    val dockerRun = DockerValidation.optional.extract(validatedRuntimeAttributes) map buildDockerRunCommand getOrElse ""
    new SharedFileSystemCommand(Seq("/bin/bash", "-c", s"cat ${jobPaths.script} | $dockerRun /bin/bash <&0"))
  }

  override lazy val stdoutSubmit = jobPaths.stdout.untailed

  override lazy val stderrSubmit = jobPaths.stderr.untailed

  override def getJob(process: Process, stdoutWriter: PathWriter, stderrWriter: PathWriter) = {
    val p = ShadowLocalBackendLifecycleActorFactory.simpleProcessFieldP.get(process)
    val pid = ShadowLocalBackendLifecycleActorFactory.unixProcessFieldPid.get(p)
    ShadowLocalJob(pid.toString)
  }

  override def getJob(exitValue: Int, stdout: Path, stderr: Path) = {
    throw new NotImplementedError("Local currently does not use a wrapper script that exits.")
  }

  override def cleanup() = {
    // TODO: This is needed until we make our own script that "returns" the pid of a child process, like qsub -terse.
    stdoutSubmit.writer.flushAndClose()
    stderrSubmit.writer.flushAndClose()
  }

  override def killArgs(job: ShadowLocalJob) = {
    SharedFileSystemCommand("kill", job.jobId)
  }

  /**
    * --rm automatically deletes the container upon exit
    * -v maps the host workflow executions directory to /root/<workflow id> on the container.
    * -i makes the run interactive, required for the cat and <&0 shenanigans that follow.
    */
  private def buildDockerRunCommand(image: String): String = {
    val dockerDir = jobPaths.callDockerRoot
    s"docker run --rm -v ${jobPaths.callRoot.toAbsolutePath}:$dockerDir -i $image"
  }
}

object ShadowLocalBackendLifecycleActorFactory {
  val simpleProcessFieldP = {
    val clazz = Class.forName("scala.sys.process.ProcessImpl.SimpleProcess")
    val field = clazz.getDeclaredField("p")
    field.setAccessible(true)
    field
  }

  val unixProcessFieldPid = {
    val clazz = Class.forName("java.lang.UNIXProcess")
    val field = clazz.getDeclaredField("pid")
    field.setAccessible(true)
    field
  }
}
