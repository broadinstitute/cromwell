package cromwell.docker

import cats.effect.IO
import cromwell.docker.DockerInfoActor.{DockerInfoContext, DockerInfoResponse}
import org.http4s.client.Client

case class MockHashResponse(hashResponse: DockerInfoResponse, nb: Int)

class DockerRegistryMock(responses: MockHashResponse*) extends DockerRegistry {
  private var responsesLeft = responses.toBuffer
  // Counts the number of elements going through this "flow"
  private var _count: Int = 0

  private def nextResponse(context: DockerInfoContext): (DockerInfoResponse, DockerInfoContext) = 
      responsesLeft.headOption match {
        case Some(mockResponse) =>
          _count += 1
          if (mockResponse.nb > 1)
            responsesLeft.update(0, mockResponse.copy(nb = mockResponse.nb - 1))
          else
            responsesLeft.remove(0)
          (mockResponse.hashResponse, context)
        // When we hit the end, loop
        case None =>
          responsesLeft = responses.toBuffer
          nextResponse(context)
  }
  
  def count() = _count

  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean = true

  override def run(dockerInfoContext: DockerInfoContext)(implicit client: Client[IO]) = IO.pure(nextResponse(dockerInfoContext))

  override def config = DockerRegistryConfig.default
}
