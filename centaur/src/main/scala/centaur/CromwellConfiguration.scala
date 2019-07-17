package centaur

import java.lang.ProcessBuilder.Redirect

import better.files.File
import com.typesafe.scalalogging.StrictLogging

trait CromwellProcess extends StrictLogging {
  def logFile: String
  def start(): Unit
  def stop(): Unit
  def isAlive: Boolean
  def cromwellConfiguration: CromwellConfiguration

  protected def runProcess(command: Array[String]): Process = {
    logger.info(s"Running: ${command.mkString(" ")}")
    val processBuilder = new java.lang.ProcessBuilder()
      .command(command: _*)
      .redirectOutput(Redirect.appendTo(File(logFile).toJava))
      .redirectErrorStream(true)
    processBuilder.start()
  }

  protected def waitProcess(process: Process, destroy: Boolean = false): Unit = {
    process.getOutputStream.flush()
    if (destroy)
      process.destroy()
    process.waitFor()
    ()
  }
}

trait CromwellConfiguration {
  def createProcess: CromwellProcess
  def logFile: String
}
