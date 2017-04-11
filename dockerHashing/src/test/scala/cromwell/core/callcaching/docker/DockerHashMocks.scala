package cromwell.core.callcaching.docker

import akka.http.scaladsl.model._
import akka.stream.scaladsl.Flow
import cromwell.core.callcaching.docker.DockerHashActor.{DockerHashContext, DockerHashResponse}
import cromwell.core.callcaching.docker.registryv2.flows.HttpFlowWithRetry.ContextWithRequest

import scala.util.Try

case class MockHttpResponse(httpResponse: Try[HttpResponse], nb: Int)

class HttpMock[T](responses: MockHttpResponse*) {
  private var responsesLeft = responses.toBuffer
  
  private def nextResponse(value: (HttpRequest, ContextWithRequest[T])): (Try[HttpResponse], ContextWithRequest[T]) = value match {
    case (request, context) =>
      responsesLeft.headOption match {
        case Some(mockResponse) =>
          if (mockResponse.nb > 1)
            responsesLeft.update(0, mockResponse.copy(nb = mockResponse.nb - 1))
          else
            responsesLeft.remove(0)
          (mockResponse.httpResponse, context)
        // When we hit the end, loop
        case None =>
          responsesLeft = responses.toBuffer
          nextResponse(value)
      }
  }
  
  def httpMock() = Flow[(HttpRequest, ContextWithRequest[T])] map nextResponse
}

case class MockHashResponse(hashResponse: DockerHashResponse, nb: Int)

class DockerFlowMock(responses: MockHashResponse*) extends DockerFlow {
  private var responsesLeft = responses.toBuffer
  // Counts the number of elements going through this "flow"
  private var _count: Int = 0

  private def nextResponse(context: DockerHashContext): (DockerHashResponse, DockerHashContext) = 
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

  override def buildFlow() = Flow[DockerHashContext] map nextResponse
  override def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash): Boolean = true
}
