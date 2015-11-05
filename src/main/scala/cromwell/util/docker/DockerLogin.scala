package cromwell.util.docker

import java.nio.charset.StandardCharsets
import java.util.Base64

import com.typesafe.config.Config
import cromwell.util.google.GoogleCredentialFactory
import lenthall.config.ScalaConfig._

/** A username and password combination. */
case class DockerLogin(username: String, password: String)

trait DockerLoginProvider {
  def dockerLogin: Option[DockerLogin]
}

class DockerHubLoginProvider(dockerTokenOption: Option[String]) extends DockerLoginProvider {
  def this(config: Config) = {
    this(config.getStringOption("docker.dockerToken"))
  }

  override val dockerLogin = dockerTokenOption map { dockerToken =>
    val bytes = Base64.getDecoder.decode(dockerToken)
    val userPass = new String(bytes, StandardCharsets.ISO_8859_1)
    userPass.indexOf(':') match {
      case -1 => DockerLogin(userPass, "")
      case index => DockerLogin(userPass.substring(0, index), userPass.substring(index + 1))
    }
  }
}

class GcrLoginProvider(config: Config) extends DockerLoginProvider {
  private lazy val credentialFactory = GoogleCredentialFactory.fromAuthScheme(config)

  override def dockerLogin = credentialFactory.freshCredential.toOption map { credential =>
    DockerLogin("_token", credential.getAccessToken)
  }
}
