package cromwell.docker.registryv2.flows.aws

import cats.effect.IO
import cromwell.docker.{DockerHashResult, DockerImageIdentifier, DockerInfoActor, DockerRegistryConfig}
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import cromwell.docker.registryv2.flows.aws.EcrUtils.{EcrForbidden, EcrNotFound, EcrUnauthorized}
import org.apache.commons.codec.digest.DigestUtils
import org.http4s.{Header, Response, Status}
import org.slf4j.{Logger, LoggerFactory}

abstract class AmazonEcrAbstract(override val config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) {
  private val logger: Logger = LoggerFactory.getLogger(this.getClass)


  /**
    * Not used as getToken is overridden
    */
  override protected def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier): String = {
    logger.warn("A call has been made to authorizationServerHostName for an ECR container but ECR auth should be done via SDK. This call will be ignored")
    ""
  }

  /**
    * Not used as getToken is overridden
    */
  override protected def buildTokenRequestHeaders(dockerInfoContext: DockerInfoActor.DockerInfoContext): List[Header] = {
    logger.warn("A call has been made to buildTokenRequestHeaders for an ECR container but ECR auth should be done via SDK. This call will be ignored")
    List.empty
  }

  /**
    * Amazon ECR repositories don't have a digest header in responses so we must make it from the manifest body
    */
  override protected def getDigestFromResponse(response: Response[IO]): IO[DockerHashResult] = response match {
    case Status.Successful(r) => digestManifest(r.bodyText)
    case Status.Unauthorized(_) => IO.raiseError(new EcrUnauthorized("Not authorized to obtain a digest from this registry. Your token may be invalid"))
    case Status.NotFound(_) => IO.raiseError(new EcrNotFound)
    case Status.Forbidden(_) => IO.raiseError(new EcrForbidden)
    case failed => failed.as[String].flatMap(body => IO.raiseError(new Exception(s"Failed to get manifest: $body")))
  }

  private def digestManifest(bodyText: fs2.Stream[IO, String]): IO[DockerHashResult] = {
    logger.info("computing sha256 digest for container manifest")
    bodyText
      .compile
      .string
      .map(data => "sha256:"+DigestUtils.sha256Hex(data))
      .map(DockerHashResult.fromString)
      .flatMap(IO.fromTry)
  }
}
