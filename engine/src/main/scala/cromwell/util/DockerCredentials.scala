package cromwell.util

import com.typesafe.config.Config
import cromwell.util.ConfigUtil._

/**
 * Encapsulate docker credential information.
 */
case class DockerCredentials(account: String, token: String)

case class DockerHubConfiguration(namespace: String, v1Registry: String, v2Registry: String)

case class DockerConfiguration(dockerCredentials: Option[DockerCredentials], dockerHubConf: DockerHubConfiguration)

/**
 * Singleton encapsulating a DockerConf instance.
 */
object DockerConfiguration {
  import lenthall.config.ScalaConfig._

  private val dockerKeys = Set("account", "token")

  def build(config: Config) = {
    val dockerConf: Option[DockerCredentials] = for {
      dockerConf <- config.getConfigOption("dockerHub")
      _ = dockerConf.warnNotRecognized(dockerKeys, "dockerHub")
      account <- dockerConf.validateString("account").toOption
      token <- dockerConf.validateString("token").toOption
    } yield new DockerCredentials(account, token)

    val dockerHubConf = {
      new DockerHubConfiguration(
        namespace = config.getStringOr("docker.hub.namespace", "docker.io"),
        v1Registry = config.getStringOr("docker.hub.v1Registry", "index.docker.io"),
        v2Registry = config.getStringOr("docker.hub.v2Registry", "registry-1.docker.io")
      )
    }
    new DockerConfiguration(dockerConf, dockerHubConf)
  }
}
