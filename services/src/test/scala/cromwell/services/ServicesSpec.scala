package cromwell.services

import akka.actor.ActorSystem
import akka.testkit.TestKit
import akka.util.Timeout
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{OneInstancePerTest, BeforeAndAfterAll, Matchers, WordSpecLike}

import scala.concurrent.duration._

abstract class ServicesSpec(serviceName: String) extends TestKit(ActorSystem(s"${serviceName}ServiceActorSpec"))
  with Matchers with WordSpecLike with BeforeAndAfterAll with ScalaFutures {

  implicit val timeout = Timeout(5.seconds)
  implicit val ec = system.dispatcher
  implicit val defaultPatience = PatienceConfig(timeout = Span(5, Seconds), interval = Span(100, Millis))

  override protected def afterAll() = {
    TestKit.shutdownActorSystem(system)
  }
}
