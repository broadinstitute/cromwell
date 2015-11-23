package cromwell.util.docker

import com.typesafe.config.ConfigFactory
import lenthall.config.ScalaConfig._

/**
  * A set of endpoints belonging to a registry.
  */
object DockerRegistry {
  def apply(endpoint: String): DockerRegistry = DockerRegistry(endpoint, endpoint)

  private val config = ConfigFactory.load()

  val DockerHub = DockerRegistry(
    config.getStringOr("docker.hub.v1Registry", "registry.hub.docker.com"),
    config.getStringOr("docker.hub.v2Registry", "registry-1.docker.io"))
 }

/**
  * A set of endpoints belonging to a registry.
  *
  * @param v1Endpoint The host for contacting the V1 API.
  * @param v2Endpoint The host for contact the V2 API.
  */
case class DockerRegistry(v1Endpoint: String, v2Endpoint: String)
