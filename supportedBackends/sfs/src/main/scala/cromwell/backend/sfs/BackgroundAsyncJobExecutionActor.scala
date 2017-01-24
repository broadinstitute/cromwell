package cromwell.backend.sfs

import java.nio.file.Path

import better.files._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.path.PathFactory._

trait BackgroundAsyncJobExecutionActor extends SharedFileSystemAsyncJobExecutionActor {

  lazy val backgroundScript: File = pathPlusSuffix(jobPaths.script, "background")

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
    val stdout = pathPlusSuffix(jobPaths.stdout, "background")
    val stderr = pathPlusSuffix(jobPaths.stderr, "background")
    val argv = Seq("/bin/bash", backgroundScript)
    new ProcessRunner(argv, stdout.path, stderr.path)
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
