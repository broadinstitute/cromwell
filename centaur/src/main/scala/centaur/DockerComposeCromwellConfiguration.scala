package centaur

import centaur.CromwellManager.ManagedCromwellPort
import com.typesafe.config.Config

object DockerComposeCromwellConfiguration {
  def apply(conf: Config): CromwellConfiguration = {
    val dockerTag = conf.getString("docker-tag")
    val dockerComposeFile = conf.getString("docker-compose-file")
    val confPath = conf.getString("conf")
    val logPath = conf.getString("log")

    new DockerComposeCromwellConfiguration(dockerTag, dockerComposeFile, confPath, logPath)
  }
}

case class DockerComposeCromwellConfiguration(dockerTag: String, dockerComposeFile: String, conf: String, logFile: String) extends CromwellConfiguration {
  override def createProcess: CromwellProcess = {
    case class DockerComposeCromwellProcess(override val cromwellConfiguration: DockerComposeCromwellConfiguration) extends CromwellProcess {

      private def composeCommand(command: String): Array[String] = {
        Array(
          "/bin/bash",
          "-c",
          s"MANAGED_CROMWELL_PORT=$ManagedCromwellPort " +
            s"CROMWELL_TAG=$dockerTag " +
            s"CROMWELL_CONFIG=$conf " +
            s"docker-compose -f $dockerComposeFile $command")
      }

      private val startCommand = composeCommand("up")
      private val stopCommand = composeCommand("down -v")
      private val rmCommand = composeCommand("rm -fsv")

      private var process: Option[Process] = None

      override def start(): Unit = {
        process = Option(runProcess(startCommand))
      }

      override def stop(): Unit = {
        waitProcess(runProcess(stopCommand))
        waitProcess(runProcess(rmCommand))
        process foreach {
          waitProcess(_, destroy = true)
        }
        process = None
      }

      override def isAlive: Boolean = process.exists {
        _.isAlive
      }

      override def logFile: String = cromwellConfiguration.logFile
    }

    DockerComposeCromwellProcess(this)
  }
}
