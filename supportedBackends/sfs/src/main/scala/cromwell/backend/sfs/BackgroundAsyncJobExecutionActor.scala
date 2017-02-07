package cromwell.backend.sfs

import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.path.Path

trait BackgroundAsyncJobExecutionActor extends SharedFileSystemAsyncJobExecutionActor {

  lazy val backgroundScript = jobPaths.script.plusExt("background")

  override def writeScriptContents(): Unit = {
    super.writeScriptContents()
    writeBackgroundScriptContents()
  }

  /**
    * Run the command via bash in the background, and echo the PID.
    */
  private def writeBackgroundScriptContents(): Unit = {
    val backgroundCommand = redirectOutputs(processArgs.argv.mkString("'", "' '", "'"))
    // $! contains the previous background command's process id (PID)
    backgroundScript.write(
      s"""|#!/bin/bash
          |BACKGROUND_COMMAND &
          |echo $$!
          |""".stripMargin.replace("BACKGROUND_COMMAND", backgroundCommand))
    ()
  }

  override def makeProcessRunner(): ProcessRunner = {
    val stdout = jobPaths.stdout.plusExt("background")
    val stderr = jobPaths.stderr.plusExt("background")
    val argv = Seq("/bin/bash", backgroundScript)
    new ProcessRunner(argv, stdout, stderr)
  }

  override def getJob(exitValue: Int, stdout: Path, stderr: Path): StandardAsyncJob = {
    val pid = stdout.contentAsString.stripLineEnd
    StandardAsyncJob(pid)
  }

  override def checkAliveArgs(job: StandardAsyncJob): SharedFileSystemCommand = {
    SharedFileSystemCommand("ps", job.jobId)
  }

  override def killArgs(job: StandardAsyncJob): SharedFileSystemCommand = {
    val killScript = jobPaths.script.plusExt("kill")
    writeKillScript(killScript, job)
    SharedFileSystemCommand("/bin/bash", killScript)
  }

  private def writeKillScript(killScript: Path, job: StandardAsyncJob): Unit = {
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
