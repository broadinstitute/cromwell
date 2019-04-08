package centaur

import java.lang.ProcessBuilder.Redirect

import better.files.File
import com.typesafe.config.Config


case class GkeCromwellConfiguration(startCommand: Array[String], stopCommand: Array[String], logFile: String) extends CromwellConfiguration {
  // Assumptions: A Cloud SQL instance has been built, a Docker image has been pushed to GCR, a GKE cluster has been built
  // and hosts 3 Cromwell deployments:
  // - a summarizer deployment with a single instance
  // - some number of frontend instances fronted by a LoadBalancer service.
  // - zero backend instances (the deployment is defined though).
  //
  // For now `start`ing and `stop`ing will be defined as scaling the backend deployment from 0 to N (where N > 0) or vice-versa.
  // A GKE configuration is quite a bit more complicated than jar or docker-compose configurations and other definitions of
  // `start` and `stop` are certainly possible, but this is what it will be for now.
  override def createProcess: CromwellProcess = {
    case class GkeCromwellProcess(override val cromwellConfiguration: GkeCromwellConfiguration) extends CromwellProcess {

      override def displayString: String = startCommand.mkString(" ")

      override def start(): Unit = {
        val processBuilder = new java.lang.ProcessBuilder()
          .command(startCommand: _*)
          .redirectOutput(Redirect.appendTo(File(logFile).toJava))
          .redirectErrorStream(true)
        processBuilder.start().waitFor()
        ()
      }

      override def stop(): Unit = {
        val processBuilder = new java.lang.ProcessBuilder()
          .command(stopCommand: _*)
          .redirectOutput(Redirect.appendTo(File(logFile).toJava))
          .redirectErrorStream(true)
        processBuilder.start().waitFor()
        ()
      }

      override def isAlive: Boolean = true

      override def logFile: String = cromwellConfiguration.logFile
    }

    GkeCromwellProcess(this)
  }
}


object GkeCromwellConfiguration {
  def apply(conf: Config): CromwellConfiguration = {
    val startCommand = conf.getString("start-command")
    val stopCommand = conf.getString("stop-command")
    val logPath = conf.getString("log")

    new GkeCromwellConfiguration(
      startCommand = Array("/bin/bash", "-c", startCommand),
      stopCommand = Array("/bin/bash", "-c", stopCommand),
      logFile = logPath
    )
  }
}
