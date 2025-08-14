package cromwell.docker

import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.docker.registryv2.flows.dockerhub.DockerHub
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.ValueReader

import scala.util.{Failure, Success}

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

/*
 * Configures a Docker mirror, which can be used to enable Cromwell to use pull-through caches for Docker images.
 * This is most useful for creating a pull-through Dockerhub cache in a cloud vendor's container registry, for example
 * Google's Artifact Registry, AWS's ECR, etc.
 *
 *  - address: The address of the mirror. The host of user-provided images will be replaced with this address.
 */
case class DockerHubMirror(address: String) extends DockerMirror with LazyLogging {
  val mirrorImage: PartialFunction[DockerImageIdentifier, DockerImageIdentifier] = {
    case i if DockerHub.isValidDockerHubHost(i.host) =>
      // if address contains slashes, terms after them should be made part of the repository, not the host name
      // Example:
      //   address: my.mirror.io/my-namespace
      //   image: docker.io/broadinstitute/cromwell
      //   mirrored image: my.mirror.io/my-namespace/broadinstitute/cromwell:latest
      //   mirrored image host: my.mirror.io
      //   mirrored image repository: my-namespace/broadinstitute
      // we can be lazy and get this behavior for free by swapping the host to the mirror address, then re-parsing the image identifier
      val mirroredImage = i.withDefaultRepository.swapHost(address)
      DockerImageIdentifier.fromString(mirroredImage.fullName) match {
        case Success(mirroredId) => mirroredId
        case Failure(exception) =>
          logger.error(
            s"Failed to mirror Docker image '${i.fullName}' using mirror address '$address'. " +
              s"Falling back to original image. Exception: ${exception.getMessage}"
          )
          i
      }
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
