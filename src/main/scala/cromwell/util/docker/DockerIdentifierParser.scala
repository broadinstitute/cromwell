package cromwell.util.docker

import scala.annotation.tailrec

/**
  * Parses a String into a DockerIdentifier.
  */
object DockerIdentifierParser {
  /**
    * Parses a String into a DockerIdentifier.
    *
    * Currently only images from Docker Hub and gcr.io are supported.
    *
    * @param identifier The identifier string.
    * @return The DockerIdentifier.
    */
  def parse(identifier: String): DockerIdentifier = parse(identifier, parsers)

  /** The list of supported parsers. */
  private val parsers = Seq(GcrDockerParser, DockerHubParser, DockerHubLibraryParser)

  /**
    * Tries the defined parsers in order.
    * If all else fails assumes the string is just a Docker Hub library, with the "latest" tag.
    *
    * @param identifier The identifier string.
    * @param parsers Remaining parsers to try.
    * @return The DockerIdentifier.
    */
  @tailrec
  private def parse(identifier: String, parsers: Seq[DockerIdentifierParser]): DockerIdentifier = {
    parsers.headOption match {
      case Some(parser) =>
        val tryMatch = parser.matchIdentifier
        if (tryMatch.isDefinedAt(identifier)) tryMatch.apply(identifier) else parse(identifier, parsers.tail)
      case None => DockerHubLibraryParser.defaultIdentifier(identifier)
    }
  }
}

/**
  * Parses a String into a DockerIdentifier.
  */
sealed trait DockerIdentifierParser {
  /**
    * @return A function that can possibly parse a String into a DockerIdentifier.
    */
  def matchIdentifier: PartialFunction[String, DockerIdentifier]
}

/**
  * Parses Docker Hub image strings that include a team/user.
  */
object DockerHubParser extends DockerIdentifierParser {
  private val hubLatest = """([^/]+)/(.*)""".r
  private val hubWithTag = """([^/]+)/(.*):([^:]+)""".r
  private val hubWithDigest = """([^/]+)/(.*)@([^@]+)""".r

  override def matchIdentifier = {
    case hubWithDigest(name, image, digest) =>
      DockerDigestIdentifier(s"$name/$image", digest, DockerRegistry.DockerHub)
    case hubWithTag(name, image, tag) =>
      DockerTagIdentifier(s"$name/$image", tag, DockerRegistry.DockerHub)
    case hubLatest(name, image) =>
      DockerTagIdentifier(s"$name/$image", "latest", DockerRegistry.DockerHub)
  }
}

/**
  * Parses Docker Hub image strings for a library, such as "ubuntu".
  */
object DockerHubLibraryParser extends DockerIdentifierParser {

  /** When all else fails, just look for the image at that name on DockerHub. */
  def defaultIdentifier(image: String) = DockerTagIdentifier(s"library/$image", "latest", DockerRegistry.DockerHub)

  private val hubLibraryWithTag = """(.*):([^:]+)""".r
  private val hubLibraryWithDigest = """(.*)@([^@]+)""".r

  override def matchIdentifier = {
    case hubLibraryWithDigest(image, digest) =>
      DockerDigestIdentifier(s"library/$image", digest, DockerRegistry.DockerHub)
    case hubLibraryWithTag(image, tag) =>
      DockerTagIdentifier(s"library/$image", tag, DockerRegistry.DockerHub)
  }
}

/**
  * Parses an image hosted on *.gcr.io.
  *
  * TODO: Shouldn't we be able to handle any host name?
  */
object GcrDockerParser extends DockerIdentifierParser {
  private val gcrLatest = """(.*\.?gcr.io)/(.*)/(.*)""".r
  private val gcrWithTag = """(.*\.?gcr.io)/(.*)/(.*):([^:]+)""".r
  private val gcrWithDigest = """(.*\.?gcr.io)/(.*)/(.*)@([^@]+)""".r

  override def matchIdentifier = {
    case gcrWithDigest(host, project, image, digest) =>
      DockerDigestIdentifier(s"$project/$image", digest, DockerRegistry(host, DockerLogin.GcrLoginOption))
    case gcrWithTag(host, project, image, tag) =>
      DockerTagIdentifier(s"$project/$image", tag, DockerRegistry(host, DockerLogin.GcrLoginOption))
    case gcrLatest(host, project, image) =>
      DockerTagIdentifier(s"$project/$image", "latest", DockerRegistry(host, DockerLogin.GcrLoginOption))
  }
}
