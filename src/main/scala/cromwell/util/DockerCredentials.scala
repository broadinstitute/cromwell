package cromwell.util

import com.typesafe.config.ConfigFactory
import cromwell.util.ConfigUtil._

/**
 * Encapsulate docker credential information.
 */
 case class DockerCredentials(account: String, token: String)

/**
 * Trait for DockerConfiguration
 */
trait DockerConfiguration {
  val dockerKeys = Set("dockerAccount", "dockerToken")

  def dockerConf: Option[DockerCredentials]
}

/**
 * Singleton encapsulating a DockerConf instance.
 */
object ProductionDockerConfiguration extends DockerConfiguration {
  override lazy val dockerConf: Option[DockerCredentials] = for {
    dockerConf <- ConfigFactory.load.getConfigOption("docker")
    _ = dockerConf.warnNotRecognized(dockerKeys, "docker")
    account <- dockerConf.validateString("dockerAccount").toOption
    token <- dockerConf.validateString("dockerToken").toOption
  } yield new DockerCredentials(account, token)
}
