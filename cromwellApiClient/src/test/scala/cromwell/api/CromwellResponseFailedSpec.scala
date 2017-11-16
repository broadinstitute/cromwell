package cromwell.api

import java.net.URL

import akka.actor.ActorSystem
import akka.http.scaladsl.model._
import akka.stream.ActorMaterializer
import akka.testkit.TestKit
import cromwell.api.model.CromwellFailedResponseException
import org.scalatest.{AsyncFlatSpecLike, BeforeAndAfterAll, Matchers}

import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.language.postfixOps

class CromwellResponseFailedSpec extends TestKit(ActorSystem()) with AsyncFlatSpecLike with Matchers with BeforeAndAfterAll {
  override def afterAll(): Unit = {
    Await.ready(system.terminate(), 1 second)
    super.afterAll()
  }
  
  implicit val materializer = ActorMaterializer()
  
  "CromwellAPIClient" should "try to fail the Future with a CromwellFailedResponseException if the HttpResponse is unsuccessful" in {
    val client = new CromwellClient(new URL("http://fakeurl"), "v1") {
      override def executeRequest(request: HttpRequest): Future[HttpResponse] = Future.successful(
        new HttpResponse(StatusCodes.ServiceUnavailable, List.empty[HttpHeader], HttpEntity(ContentTypes.`application/json`,
          """{
            |  "status": "fail",
            |  "message": "Cromwell service shutting down"
            |}
          """.stripMargin), HttpProtocols.`HTTP/1.1`)
      )
    }

    recoverToExceptionIf[CromwellFailedResponseException] { client.version(scala.concurrent.ExecutionContext.global) } map { exception =>
      assert(exception.status == "fail")
      assert(exception.message == "Cromwell service shutting down")
    }
  }
}
