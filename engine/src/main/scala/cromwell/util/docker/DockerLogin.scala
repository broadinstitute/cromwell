package cromwell.util.docker

import java.nio.charset.StandardCharsets
import java.util.Base64

import com.google.api.client.auth.oauth2.Credential
import com.typesafe.config.Config
import cromwell.filesystems.gcs.GoogleAuthMode
import GoogleAuthMode._
import lenthall.config.ScalaConfig._

/** A username and password combination. */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class DockerLogin(username: String, password: String)

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
trait DockerLoginProvider {
  def dockerLogin: Option[DockerLogin]
}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
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

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class GcrLoginProvider(credential: Credential) extends DockerLoginProvider {

  override def dockerLogin = credential.freshCredential.toOption map { credential =>
    DockerLogin("_token", credential.getAccessToken)
  }
}
