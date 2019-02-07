package cromwell.docker

import cats.effect.IO
import cromwell.docker.DockerInfoActor.{DockerInfoContext, DockerInfoResponse}
import org.http4s.client.Client

/**
  * Interface used by the docker hash actor to build a flow and validate whether or not it can accept an image.
  */
trait DockerRegistry {
  def run(dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]): IO[(DockerInfoResponse, DockerInfoContext)]
  def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean
  def config: DockerRegistryConfig
}
