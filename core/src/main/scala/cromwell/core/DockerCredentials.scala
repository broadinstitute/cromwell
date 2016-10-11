package cromwell.core

import com.typesafe.config.Config
import cromwell.core.ConfigUtil._

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

  private val dockerKeys = Set("account", "token")

  def build(config: Config) = {
    import net.ceedubs.ficus.Ficus._
    val dockerConf: Option[DockerCredentials] = for {
      dockerConf <- config.as[Option[Config]]("dockerhub")
      _ = dockerConf.warnNotRecognized(dockerKeys, "dockerhub")
      account <- dockerConf.validateString("account").toOption
      token <- dockerConf.validateString("token").toOption
    } yield DockerCredentials(account, token)

    val dockerHubConf = {
      DockerHubConfiguration(
        namespace = config.as[Option[String]]("docker.hub.namespace").getOrElse("docker.io"),
        v1Registry = config.as[Option[String]]("docker.hub.v1Registry").getOrElse("index.docker.io"),
        v2Registry = config.as[Option[String]]("docker.hub.v2Registry").getOrElse("registry-1.docker.io")
      )
    }
    new DockerConfiguration(dockerConf, dockerHubConf)
  }
}
