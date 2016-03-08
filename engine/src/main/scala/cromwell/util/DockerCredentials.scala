package cromwell.util

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.util.ConfigUtil._

import scala.util.Try

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

  private val dockerKeys = Set("dockerAccount", "dockerToken")

  lazy val dockerConf = build(ConfigFactory.load)

  def build(conf: Config) = {
    val dockerConf: Option[DockerCredentials] = for {
      dockerConf <- conf.getConfigOption("docker")
      _ = dockerConf.warnNotRecognized(dockerKeys, "docker")
      account <- dockerConf.validateString("dockerAccount").toOption
      token <- dockerConf.validateString("dockerToken").toOption
    } yield new DockerCredentials(account, token)

    val dockerHubConf = {
      new DockerHubConfiguration(
        namespace = conf.getStringOr("docker.hub.namespace", "docker.io"),
        v1Registry = conf.getStringOr("docker.hub.v1Registry", "index.docker.io"),
        v2Registry = conf.getStringOr("docker.hub.v2Registry", "registry-1.docker.io")
      )
    }
    new DockerConfiguration(dockerConf, dockerHubConf)
  }

}
