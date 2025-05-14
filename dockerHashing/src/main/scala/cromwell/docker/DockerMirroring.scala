package cromwell.docker

import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.docker.registryv2.flows.dockerhub.DockerHub
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.ValueReader

case class DockerMirroring(mirrors: List[DockerMirror]) {
  def mirrorImage(image: DockerImageIdentifier): Option[DockerImageIdentifier] =
    mirrors.flatMap(_.mirrorImage.lift(image)).headOption
}

object DockerMirroring {
  def fromConfig(config: Config): Option[DockerMirroring] = {
    val mirrors = config
      .getAs[Config]("docker-mirror")
      .map { mirrorConfig =>
        val dockerhubMirror = mirrorConfig.getAs[DockerHubMirror]("dockerhub")
        // Add support for additional repositories here

        dockerhubMirror.toList
      }
      .getOrElse(List[DockerMirror]())
    Option.when(mirrors.nonEmpty)(DockerMirroring(mirrors))
  }
}

sealed trait DockerMirror {
  val mirrorImage: PartialFunction[DockerImageIdentifier, DockerImageIdentifier]
}

case class DockerHubMirror(address: String) extends DockerMirror {
  val mirrorImage: PartialFunction[DockerImageIdentifier, DockerImageIdentifier] = {
    case i if DockerHub.isValidDockerHubHost(i.host) => i.swapHost(address)
  }
}

object DockerHubMirror extends LazyLogging {
  implicit val dockerHubMirrorOptionValueReader: ValueReader[Option[DockerHubMirror]] =
    (config: Config, path: String) =>
      config.getAs[Config](path) flatMap { dockerMirrorConfig =>
        val enabled = dockerMirrorConfig.as[Boolean]("enabled")
        val address = dockerMirrorConfig.as[Option[String]]("address").getOrElse("")
        if (address.isEmpty)
          logger.warn(
            "Potential misconfiguration: docker-mirror.dockerhub.enabled=true with no address provided. " +
              "Mirroring will be disabled."
          )
        Option.when(enabled && address.nonEmpty)(DockerHubMirror(address))
      }
}
