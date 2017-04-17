package cromwell.util.docker

import com.typesafe.config.ConfigFactory
import lenthall.config.ScalaConfig._

/**
  * A set of endpoints belonging to a registry.
  */
object DockerRegistry {
  def apply(host: String, login: Option[DockerLogin]): DockerRegistry = DockerRegistry(host, host, host, login)

  private val config = ConfigFactory.load()

  val DockerHub = DockerRegistry(
    config.getStringOr("docker.hub.namespace", "docker.io"),
    config.getStringOr("docker.hub.v1Registry", "registry.hub.docker.com"),
    config.getStringOr("docker.hub.v2Registry", "registry-1.docker.io"),
    None)
 }

/**
  * A set of endpoints belonging to a registry.
  *
  * @param v1Hostname The host for contacting the V1 API.
  * @param v2Hostname The host for contacting the V2 API.
  */
case class DockerRegistry(namespace: String, v1Hostname: String, v2Hostname: String, login: Option[DockerLogin])
