package cromwell.util.docker

import com.typesafe.config.{Config, ConfigFactory}

import scala.annotation.tailrec

/**
  * Parses a String into a DockerIdentifier.
  *
  * @param config The configuration for the parsers.
  */
class DockerIdentifierParser(config: Config) {
  private lazy val dockerHubRegistryParser = new DockerHubRegistryParser(config)
  private lazy val googleContainerRegistryParser = new GoogleContainerRegistryParser(config)

  /** The list of supported registry parsers. */
  private lazy val registryParsers = Seq(googleContainerRegistryParser, dockerHubRegistryParser)

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
  lazy val Default = new DockerIdentifierParser(ConfigFactory.load)
}
