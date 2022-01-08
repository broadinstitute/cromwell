package centaur

import java.lang.ProcessBuilder.Redirect

import better.files.File
import com.typesafe.scalalogging.StrictLogging
import scala.collection.JavaConverters._

trait CromwellProcess extends StrictLogging {
  def logFile: String
  def start(): Unit
  def stop(): Unit
  def isAlive: Boolean
  def cromwellConfiguration: CromwellConfiguration

  protected def runProcess(command: Array[String], additionalEnv: Map[String, String]): Process = {
    logger.info(s"Running: ${command.mkString(" ")}")
    val processBuilder = new java.lang.ProcessBuilder()
      .command(command: _*)
      .redirectOutput(Redirect.appendTo(File(logFile).toJava))
      .redirectErrorStream(true)
    processBuilder.environment().putAll(additionalEnv.asJava)
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
