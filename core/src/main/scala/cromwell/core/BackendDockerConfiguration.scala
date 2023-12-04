package cromwell.core

import java.util.Base64

import com.typesafe.config.Config
import cromwell.core.ConfigUtil._

import scala.util.Try

object DockerCredentials {
  def unapply(arg: DockerCredentials): Option[String] = Option(arg.token)
}

object DockerCredentialUsernameAndPassword {
  private val tokenStringFormat = raw"([^:]*):(.*)".r

  def unapply(arg: DockerCredentials): Option[(String, String)] = Try(
    new String(Base64.getDecoder.decode(arg.token))
  ).toOption match {
    case Some(tokenStringFormat(username, password)) => Some((username, password))
    case _ => None
  }
}

/**
  * Encapsulate docker credential information.
  */
class DockerCredentials(val token: String, val keyName: Option[String], val authName: Option[String])

case class BackendDockerConfiguration(dockerCredentials: Option[DockerCredentials])

/**
  * Singleton encapsulating a DockerConf instance.
  */
object BackendDockerConfiguration {

  private val dockerKeys = Set("token", "auth", "key-name")

  def build(config: Config) = {
    import net.ceedubs.ficus.Ficus._
    val dockerConf: Option[DockerCredentials] = for {
      dockerConf <- config.as[Option[Config]]("dockerhub")
      _ = dockerConf.warnNotRecognized(dockerKeys, "dockerhub")
      token <- dockerConf.validateString("token").toOption
      authName = dockerConf.as[Option[String]]("auth")
      keyName = dockerConf.as[Option[String]]("key-name")
    } yield new DockerCredentials(token = token, authName = authName, keyName = keyName)

    new BackendDockerConfiguration(dockerConf)
  }
}
