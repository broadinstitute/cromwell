package cromwell.backend.sfs

import java.nio.file.Path

import better.files._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.path.PathFactory._

trait BackgroundAsyncJobExecutionActor extends SharedFileSystemAsyncJobExecutionActor {

  override def makeProcessRunner(): ProcessRunner = {
    val backgroundScript = pathPlusSuffix(jobPaths.script, "background")
    writeBackgroundScript(backgroundScript, processArgs.argv.mkString("'", "' '", "'"))
    val stdout = pathPlusSuffix(jobPaths.stdout, "background")
    val stderr = pathPlusSuffix(jobPaths.stderr, "background")
    val argv = Seq("/bin/bash", backgroundScript)
    new ProcessRunner(argv, stdout.path, stderr.path)
  }

  private def writeBackgroundScript(backgroundScript: File, backgroundCommand: String): Unit = {
    /*
    Run the `backgroundCommand` in the background. Redirect the stdout and stderr to the appropriate files. While not
    necessary, mark the job as not receiving any stdin by pointing it at /dev/null.

    If the `backgroundCommand` errors for some reason, put a "-1" into the rc file.

    Finally, run all of the above in the bash background, and return the PID of the backgrounded command.

    bashism | english
    --------|--------------------------------------------------------------------------
       >    | redirect stdout to <file>
       2>   | redirect stderr to <file>
       <    | redirect stdin from <file>
       ||   | if the previous command fails, then run the following command
       >    | redirect stdout to <file>
       &    | send the entire compound command, including the || to the background
       $!   | a variable containing the previous background command's process id (PID)
     */
    backgroundScript.write(
      s"""|#!/bin/bash
          |$backgroundCommand \\
          |  > ${jobPaths.stdout} \\
          |  2> ${jobPaths.stderr} \\
          |  < /dev/null \\
          |  || echo -1 \\
          |  > ${jobPaths.returnCode} \\
          |  &
          |echo $$!
          |""".stripMargin)
    ()
  }

  override def getJob(exitValue: Int, stdout: Path, stderr: Path): StandardAsyncJob = {
    val pid = File(stdout).contentAsString.stripLineEnd
    StandardAsyncJob(pid)
  }

  override def checkAliveArgs(job: StandardAsyncJob): SharedFileSystemCommand = {
    SharedFileSystemCommand("ps", job.jobId)
  }

  override def killArgs(job: StandardAsyncJob): SharedFileSystemCommand = {
    val killScript = pathPlusSuffix(jobPaths.script, "kill")
    writeKillScript(killScript, job)
    SharedFileSystemCommand("/bin/bash", killScript)
  }

  private def writeKillScript(killScript: File, job: StandardAsyncJob): Unit = {
    /*
    Use pgrep to find the children of a process, and recursively kill the children before killing the parent.
     */
    killScript.write(
      s"""|#!/bin/bash
          |kill_children() {
          |  local pid=$$1
          |  for cpid in $$(pgrep -P $$pid); do
          |    kill_children $$cpid
          |  done
          |  echo killing $$pid
          |  kill $$pid
          |}
          |
          |kill_children ${job.jobId}
          |""".stripMargin)
    ()
  }
}
