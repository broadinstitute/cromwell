package cromwell.azure_e2e

import akka.actor.ActorSystem
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.BeforeAndAfterAll
import org.scalatest.matchers.should.Matchers._
import akka.http.scaladsl.Http
import akka.http.scaladsl.model._
import org.scalatest.matchers.should.Matchers

import scala.concurrent.Future
import scala.util.{Failure, Success}

class AzureRunWorkflowSpec extends AnyFlatSpec with BeforeAndAfterAll with Matchers {

//  implicit val system = ActorSystem()

  //Before building out the test, confirm that env variables defined in the GHA is being passed into this testing env
  "Run workflow" should "be successfully submitted with Cromwell on Azure" in {
      val owner: String = sys.env("OWNER")
      owner should be ("hermione.owner@quality.firecloud.org")

//    val responseFuture: Future[HttpResponse] = Http().singleRequest(
//      HttpRequest(
//
//      )
//    )
//    responseFuture onComplete {
//      case Success(res) => println("test")
//      case Failure(_) => println("something wrong")
//    }
  }
}
