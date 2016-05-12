package cromwell.util.docker

import com.google.api.client.auth.oauth2.Credential
import cromwell.core.DockerConfiguration

/**
  * Parses a String into a DockerIdentifier for a particular registry.
  */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
sealed trait DockerRegistryIdentifierParser {
  /**
    * @return A function that can possibly parse a String into a DockerIdentifier.
    */
  def matchIdentifier: PartialFunction[String, DockerIdentifier]
}

/**
  * Parses Docker Hub image strings.
  */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class DockerHubRegistryParser(dockerConfig: DockerConfiguration) extends DockerRegistryIdentifierParser {
  private lazy val dockerHubLoginProvider = new DockerHubLoginProvider(dockerConfig.dockerCredentials map { _.token })

  private lazy val dockerHubRegistry = DockerRegistry(
    dockerConfig.dockerHubConf.namespace,
    dockerConfig.dockerHubConf.v1Registry,
    dockerConfig.dockerHubConf.v2Registry,
    dockerHubLoginProvider
  )
  /** When all else fails, just look for the image at that name on DockerHub. */
  def defaultIdentifier(image: String) = DockerTagIdentifier(s"library/$image", "latest", dockerHubRegistry)

  private val hubLatest = """([^/]+)/(.*)""".r
  private val hubWithTag = """([^/]+)/(.*):([^:]+)""".r
  private val hubWithDigest = """([^/]+)/(.*)@([^@]+)""".r
  private val hubLibraryWithTag = """(.*):([^:]+)""".r
  private val hubLibraryWithDigest = """(.*)@([^@]+)""".r

  override def matchIdentifier = {
    case hubWithDigest(name, image, digest) => DockerDigestIdentifier(s"$name/$image", digest, dockerHubRegistry)
    case hubWithTag(name, image, tag) => DockerTagIdentifier(s"$name/$image", tag, dockerHubRegistry)
    case hubLatest(name, image) => DockerTagIdentifier(s"$name/$image", "latest", dockerHubRegistry)
    case hubLibraryWithDigest(image, digest) => DockerDigestIdentifier(s"library/$image", digest, dockerHubRegistry)
    case hubLibraryWithTag(image, tag) => DockerTagIdentifier(s"library/$image", tag, dockerHubRegistry)
  }
}

/**
  * Parses an image hosted on *.gcr.io.
  */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class GoogleContainerRegistryParser(credential: Credential) extends DockerRegistryIdentifierParser {
  private val gcrLoginProvider = new GcrLoginProvider(credential)
  private val gcrLatest = """(.*\.?gcr.io)/(.*)/(.*)""".r
  private val gcrWithTag = """(.*\.?gcr.io)/(.*)/(.*):([^:]+)""".r
  private val gcrWithDigest = """(.*\.?gcr.io)/(.*)/(.*)@([^@]+)""".r

  override def matchIdentifier = {
    case gcrWithDigest(host, project, image, digest) =>
      DockerDigestIdentifier(s"$project/$image", digest, DockerRegistry(host, gcrLoginProvider))
    case gcrWithTag(host, project, image, tag) =>
      DockerTagIdentifier(s"$project/$image", tag, DockerRegistry(host, gcrLoginProvider))
    case gcrLatest(host, project, image) =>
      DockerTagIdentifier(s"$project/$image", "latest", DockerRegistry(host, gcrLoginProvider))
  }
}
