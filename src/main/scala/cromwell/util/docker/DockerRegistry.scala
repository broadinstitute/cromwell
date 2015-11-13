package cromwell.util.docker

import com.typesafe.config.ConfigFactory
import lenthall.config.ScalaConfig._

object DockerRegistry {
  def apply(endpoint: String): DockerRegistry = DockerRegistry(endpoint, endpoint)

  private val config = ConfigFactory.load()

  val DockerHub = DockerRegistry(
    config.getStringOr("docker.hub.v1Registry", "registry.hub.docker.com"),
    config.getStringOr("docker.hub.v2Registry", "registry-1.docker.io"))
 }

case class DockerRegistry(v1Endpoint: String, v2Endpoint: String)
