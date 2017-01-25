package cromwell.server

import akka.pattern.ask
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import org.scalatest.concurrent.{PatienceConfiguration, ScalaFutures}
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import spray.http.{ContentTypes, HttpRequest, HttpResponse, Timedout}

import scala.concurrent.duration._

class CromwellServerSpec extends FlatSpec with Matchers with Mockito with ScalaFutures {
  implicit val timeout: Timeout = 5.seconds

  it should "return 500 errors as Json" in {
    val cromwellSystem = new CromwellSystem {}

    val cromwellServerActor = cromwellSystem.actorSystem.actorOf(CromwellServerActor.props(ConfigFactory.empty()))
    val response = cromwellServerActor.ask(Timedout(mock[HttpRequest]))

    response.futureValue(PatienceConfiguration.Timeout(timeout.duration)) match {
      case h: HttpResponse =>
        h.entity.toOption match {
          case Some(e) => e.contentType.toString() shouldBe ContentTypes.`application/json`.mediaType.value.toString
          case None => fail()
        }
      case _ => fail()
    }

    cromwellSystem.shutdownActorSystem()
  }
}
