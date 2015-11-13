package cromwell.util.docker

import scala.annotation.tailrec

object DockerIdentifierParser {
  def parse(identifier: String): DockerIdentifier = parse(identifier, parsers)

  private val parsers = Seq(GcrDockerParser, DockerHubParser, DockerHubLibraryParser)

  @tailrec
  private def parse(identifier: String, parsers: Seq[DockerIdentifierParser]): DockerIdentifier = {
    parsers.headOption match {
      case Some(parser) =>
        val tryMatch = parser.matchIdentifier
        if (tryMatch.isDefinedAt(identifier)) tryMatch.apply(identifier) else parse(identifier, parsers.tail)
      case None => throw new Error("Docker Hub Library Parser should have matched.")
    }
  }
}

sealed trait DockerIdentifierParser {
  def matchIdentifier: PartialFunction[String, DockerIdentifier]
}

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

object DockerHubLibraryParser extends DockerIdentifierParser {
  private val hubLibraryLatest = """.*""".r
  private val hubLibraryWithTag = """(.*):([^:]+)""".r
  private val hubLibraryWithDigest = """(.*)@([^@]+)""".r

  override def matchIdentifier = {
    case hubLibraryWithDigest(image, digest) =>
      DockerDigestIdentifier(s"library/$image", digest, DockerRegistry.DockerHub)
    case hubLibraryWithTag(image, tag) =>
      DockerTagIdentifier(s"library/$image", tag, DockerRegistry.DockerHub)
    case hubLibraryLatest(image) =>
      DockerTagIdentifier(s"library/$image", "latest", DockerRegistry.DockerHub)
  }
}

object GcrDockerParser extends DockerIdentifierParser {
  private val gcrLatest = """(.*\.?gcr.io)/(.*)/(.*)""".r
  private val gcrWithTag = """(.*\.?gcr.io)/(.*)/(.*):([^:]+)""".r
  private val gcrWithDigest = """(.*\.?gcr.io)/(.*)/(.*)@([^@]+)""".r

  override def matchIdentifier = {
    case gcrWithDigest(host, project, image, digest) =>
      DockerDigestIdentifier(s"$project/$image", digest, DockerRegistry(host))
    case gcrWithTag(host, project, image, tag) =>
      DockerTagIdentifier(s"$project/$image", tag, DockerRegistry(host))
    case gcrLatest(host, project, image) =>
      DockerTagIdentifier(s"$project/$image", "latest", DockerRegistry(host))
  }
}
