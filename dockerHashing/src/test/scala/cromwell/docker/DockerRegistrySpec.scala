package cromwell.docker

import akka.actor.{ActorRef, Scheduler}
import akka.testkit.ImplicitSender
import cats.effect.{ContextShift, IO}
import cromwell.core.TestKitSuite

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

abstract class DockerRegistrySpec extends TestKitSuite with ImplicitSender {
  implicit val executionContext: ExecutionContext = system.dispatcher
  implicit val contextShift: ContextShift[IO] = IO.contextShift(executionContext)
  implicit val scheduler: Scheduler = system.scheduler

  protected def registryFlows: Seq[DockerRegistry]

  // Disable cache by setting a cache size of 0 - A separate test tests the cache
  lazy val dockerActor: ActorRef = system.actorOf(
    props = DockerInfoActor.props(registryFlows, 1000, 20.minutes, 0),
    name = "dockerActor",
  )

  def dockerImage(string: String): DockerImageIdentifier = DockerImageIdentifier.fromString(string).get

  def makeRequest(string: String): DockerInfoRequest = {
    DockerInfoRequest(dockerImage(string))
  }

  override protected def afterAll(): Unit = {
    system.stop(dockerActor)
    super.afterAll()
  }
}
