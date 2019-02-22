package cromwell.docker

import akka.testkit.ImplicitSender
import cats.effect.IO
import cromwell.core.TestKitSuite

import scala.concurrent.duration._

abstract class DockerRegistrySpec(actorSystemName: String) extends TestKitSuite(actorSystemName) with ImplicitSender {
  implicit val ex = system.dispatcher
  implicit val cs = IO.contextShift(ex)
  implicit val scheduler = system.scheduler

  protected def registryFlows: Seq[DockerRegistry]

  // Disable cache by setting a cache size of 0 - A separate test tests the cache
  lazy val dockerActor = system.actorOf(DockerInfoActor.props(registryFlows, 1000, 20.minutes, 0))

  def dockerImage(string: String) = DockerImageIdentifier.fromString(string).get

  def makeRequest(string: String) = {
    DockerInfoRequest(dockerImage(string))
  }

  override protected def afterAll() = {
    system.stop(dockerActor)
    super.afterAll()
  }
}
