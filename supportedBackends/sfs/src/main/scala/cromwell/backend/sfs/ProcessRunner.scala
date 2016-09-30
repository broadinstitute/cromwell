package cromwell.backend.sfs

import java.nio.file.Path

/**
  * Runs a process and sends the stdout and stderr to a file path.
  *
  * @param argv The command to run plus arguments.
  * @param stdoutPath The path to the stdout.
  * @param stderrPath The path to the stderr.
  */
class ProcessRunner(val argv: Seq[Any], val stdoutPath: Path, val stderrPath: Path) {
  def run(): Int = {
    // NOTE: Scala's SimpleProcess.exitValue() calls outputThreads.foreach(_.join()) that apparently blocks until sub
    // processes finish. WE DO NOT WANT TO WAIT FOR SUB PROCESSES TO FINISH, especially when trying to background a job!
    val processBuilder = new java.lang.ProcessBuilder()
    processBuilder.command(argv.map(_.toString): _*)
    processBuilder.redirectOutput(stdoutPath.toFile)
    processBuilder.redirectError(stderrPath.toFile)
    val process = processBuilder.start()
    process.waitFor()
  }
}
