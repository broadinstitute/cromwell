package cromwell.docker

import com.typesafe.config.Config
import cromwell.docker.registryv2.flows.dockerhub.DockerHub
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.ValueReader

case class DockerMirroring(mirrors: List[DockerMirror]) {
  def mirrorImage(image: DockerImageIdentifier): Option[DockerImageIdentifier] =
    mirrors.collectFirst(_.mirrorImage(image))
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
    if (mirrors.nonEmpty)
      Option(DockerMirroring(mirrors))
    else None
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

object DockerHubMirror {
  implicit val dockerHubMirrorOptionValueReader: ValueReader[Option[DockerHubMirror]] =
    (config: Config, path: String) =>
      config.getAs[Config](path) flatMap { dockerMirrorConfig =>
        val enabled = dockerMirrorConfig.as[Boolean]("enabled")
        // TODO how worried are we about enabled=true with an empty string?
        // TODO basically equates to no mirroring
        val address = dockerMirrorConfig.as[Option[String]]("address").getOrElse("")
        if (enabled)
          Option(DockerHubMirror(address))
        else None
      }
}
