package cromwell.util.docker

import com.google.api.client.auth.oauth2.Credential
import cromwell.util.DockerConfiguration
import cromwell.util.google.GoogleCredentialFactory

import scala.annotation.tailrec
import scala.util.Try

/**
  * Parses a String into a DockerIdentifier.
  */
class DockerIdentifierParser(dockerConf: DockerConfiguration, googleCredentials: Option[Credential]) {

  private lazy val dockerHubRegistryParser = new DockerHubRegistryParser(dockerConf)
  private lazy val googleContainerRegistryParser = googleCredentials map { new GoogleContainerRegistryParser(_) }

  /** The list of supported registry parsers. */
  private lazy val registryParsers = Seq(googleContainerRegistryParser, Option(dockerHubRegistryParser)).flatten

  /**
    * Parses a String into a DockerIdentifier.
    *
    * Currently only images from Docker Hub and gcr.io are supported.
    *
    * @param identifier The identifier string.
    * @return The DockerIdentifier.
    */
  def parse(identifier: String): DockerIdentifier = parse(identifier, registryParsers)

  /**
    * Tries the defined parsers in order.
    * If all else fails assumes the string is just a Docker Hub library, with the "latest" tag.
    *
    * @param identifier The identifier string.
    * @param parsers Remaining registries to try.
    * @return The DockerIdentifier.
    */
  @tailrec
  private def parse(identifier: String, parsers: Seq[DockerRegistryIdentifierParser]): DockerIdentifier = {
    parsers.headOption match {
      case Some(parser) =>
        val tryMatch = parser.matchIdentifier
        if (tryMatch.isDefinedAt(identifier)) tryMatch.apply(identifier) else parse(identifier, parsers.tail)
      case None => dockerHubRegistryParser.defaultIdentifier(identifier)
    }
  }
}

object DockerIdentifierParser {
  lazy val Default = new DockerIdentifierParser(DockerConfiguration.dockerConf, Try(GoogleCredentialFactory.fromCromwellAuthScheme).toOption)
}
