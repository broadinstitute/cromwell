package cromwell.core.callcaching.docker

import akka.http.scaladsl.Http
import akka.stream.ActorMaterializer
import akka.testkit.ImplicitSender
import cromwell.core.TestKitSuite
import cromwell.core.callcaching.docker.DockerHashActor._
import cromwell.core.callcaching.docker.registryv2.flows.HttpFlowWithRetry.ContextWithRequest

import scala.concurrent.duration._

abstract class DockerFlowSpec(actorSystemName: String) extends TestKitSuite(actorSystemName) with ImplicitSender {
  implicit val materializer = ActorMaterializer()
  implicit val ex = system.dispatcher
  implicit val scheduler = system.scheduler

  lazy val httpPool = Http().superPool[ContextWithRequest[DockerHashContext]]()

  protected def registryFlows: Seq[DockerFlow]

  // Disable cache by setting a cache size of 0 - A separate test tests the cache
  lazy val dockerActor = system.actorOf(DockerHashActor.props(registryFlows, 1000, 20.minutes, 0)(materializer))

  def dockerImage(string: String) = DockerImageIdentifier.fromString(string).get.asInstanceOf[DockerImageIdentifierWithoutHash]

  def makeRequest(string: String) = {
    DockerHashRequest(dockerImage(string))
  }

  override protected def afterAll() = {
    system.stop(dockerActor)
    Http().shutdownAllConnectionPools()
    materializer.shutdown()
    super.afterAll()
  }
}
