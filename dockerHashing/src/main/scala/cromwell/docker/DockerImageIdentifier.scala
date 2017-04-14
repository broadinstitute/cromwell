package cromwell.docker

import scala.util.{Failure, Success, Try}

sealed trait DockerImageIdentifier {
  def host: Option[String]
  def repository: String
  def image: String
  def reference: String

  lazy val name = s"$repository/$image"
  lazy val hostAsString = host map { h => s"$h/" } getOrElse ""
  lazy val fullName = s"$hostAsString$repository/$image:$reference"
}

case class DockerImageIdentifierWithoutHash(host: Option[String], repository: String, image: String, reference: String) extends DockerImageIdentifier {
  def withHash(hash: DockerHashResult) =  DockerImageIdentifierWithHash(host, repository, image, reference, hash)
}

case class DockerImageIdentifierWithHash(host: Option[String], repository: String, image: String, reference: String, hash: DockerHashResult) extends DockerImageIdentifier {
  override lazy val fullName: String = s"$hostAsString$repository/$image@${hash.algorithmAndHash}"
}

object DockerImageIdentifier {
  // See https://github.com/docker-library/official-images/tree/master/library
  private val DefaultDockerRepo = "library"
  private val DefaultDockerTag = "latest"
  private val TagSeparator = ":"
  private val DigestSeparator = "@"
  
  private val DockerStringRegex =
    s"""
       (?x)                                     # Turn on comments and whitespace insensitivity
       
       (                                        # Begin capturing group for name
          [a-z0-9]+(?:[._-][a-z0-9]+)*          # API v2 name component regex - see https://docs.docker.com/registry/spec/api/#/overview
          (?:/[a-z0-9]+(?:[._-][a-z0-9]+)*)*    # Optional additional name components separated by /
       )                                        # End capturing group for name
       
       (?:   
          (                                     # Begin capturing group for tag separator
            :|@                                 # Tag separator. ':' is followed by a tag, '@' is followed by a digest
          )                                     # End capturing group for tag separator
       
          (                                     # Begin capturing group for reference 
            [A-Za-z0-9]+(?:[.:_-][A-Za-z0-9]+)* # Reference
          )                                     # End capturing group for reference  
       )?                                      
       """.trim.r
  
  def fromString(dockerString: String): Try[DockerImageIdentifier] = {
    dockerString match {
        // Just a name means latest tag implicitly
      case DockerStringRegex(name, null, null) => buildId(name, TagSeparator,  DefaultDockerTag)
      case DockerStringRegex(name, tagSeparator, reference) => buildId(name, tagSeparator, reference)
      case _ => Failure(new IllegalArgumentException(s"Docker image $dockerString has an invalid syntax."))
    }
  }
  
  private def isRegistryHostName(str: String) = str.contains('.')
  
  private def buildId(name: String, tagSeparator: String, reference: String) = {
    val (dockerHost, dockerRepo, dockerImage) = name.split('/').toList match {
      // If just one component (e.g ubuntu), assume default repo
      case image :: Nil => (None, DefaultDockerRepo, image)
      // If repo/image (e.g broadinstitute/cromwell) without host
      case repo :: image :: Nil if !isRegistryHostName(repo) => (None, repo, image)
      // If host/image (e.g index.docker.io/ubuntu), assume default repo
      case host :: image :: Nil if isRegistryHostName(host) => (Option(host), DefaultDockerRepo, image)
      // Not a host followed more than one components
      case nothost :: rest if !isRegistryHostName(nothost) => (None, s"$nothost/${rest.init.mkString("/")}", rest.last)
      // A host followed more than one components (e.g gcr.io/google-containers/alpine-with-bash)
      case host :: rest if isRegistryHostName(host) => (Option(host), rest.init.mkString("/"), rest.last)
    }
    
    tagSeparator match {
      case DigestSeparator => Try(DockerHashResult(reference)) map { hash => DockerImageIdentifierWithHash(dockerHost, dockerRepo, dockerImage, reference, hash) }
      case TagSeparator => Success(DockerImageIdentifierWithoutHash(dockerHost, dockerRepo, dockerImage, reference))
        // Should have been caught by the regex, but in case..
      case other => Failure(new IllegalArgumentException(s"Invalid separator between image and tag in docker attribute: $other"))
    }
  }
}
