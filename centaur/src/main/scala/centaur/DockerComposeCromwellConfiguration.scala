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

      private def composeCommand(command: String*): Array[String] = {
        Array("docker-compose", "-f", dockerComposeFile) ++ command
      }

      private val startCommand = composeCommand("up", "--abort-on-container-exit")
      private val logsCommand = composeCommand("logs")
      private val stopCommand = composeCommand("down", "-v")
      private val rmCommand = composeCommand("rm", "-fsv")
      private val envVariables = Map[String, String](
        "CROMWELL_BUILD_CENTAUR_MANAGED_PORT" -> ManagedCromwellPort.toString,
        "CROMWELL_BUILD_CENTAUR_MANAGED_TAG" -> dockerTag,
        "CROMWELL_BUILD_CENTAUR_MANAGED_CONFIG" -> conf,
      )

      private var process: Option[Process] = None

      override def start(): Unit = {
        process = Option(runProcess(startCommand, envVariables))
      }

      override def stop(): Unit = {
        if (!isAlive) {
          // When `docker-compose up` starts successfully it attaches to the containers and prints the logs. But when
          // `docker-compose up` fails to start and exits prematurely then we need to manually retrieve the logs.
          waitProcess(runProcess(logsCommand, envVariables))
        }
        waitProcess(runProcess(stopCommand, envVariables))
        waitProcess(runProcess(rmCommand, envVariables))
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
