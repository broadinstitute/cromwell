package cromwell.core.callcaching.docker

import akka.NotUsed
import akka.stream.{FlowShape, Graph}

/**
  * Interface used by the docker hash actor to build a flow and validate whether or not it can accept an image.
  */
trait DockerFlow {
  def buildFlow(): Graph[FlowShape[DockerHashActor.DockerHashContext, (DockerHashActor.DockerHashResponse, DockerHashActor.DockerHashContext)], NotUsed]
  def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash): Boolean
}
