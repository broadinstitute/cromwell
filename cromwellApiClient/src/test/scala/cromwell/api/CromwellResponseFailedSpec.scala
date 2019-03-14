package cromwell.api

import java.net.URL

import akka.actor.ActorSystem
import akka.http.scaladsl.model._
import akka.http.scaladsl.unmarshalling.Unmarshal
import akka.stream.ActorMaterializer
import akka.testkit._
import cromwell.api.CromwellClient.UnsuccessfulRequestException
import org.scalatest.{AsyncFlatSpecLike, BeforeAndAfterAll, Matchers}

import scala.concurrent.duration._
import scala.concurrent.{Await, Future}

class CromwellResponseFailedSpec extends TestKit(ActorSystem()) with AsyncFlatSpecLike with Matchers with BeforeAndAfterAll {
  override def afterAll(): Unit = {
    Await.ready(system.terminate(), 10.seconds.dilated)
    super.afterAll()
  }
  
  implicit val materializer = ActorMaterializer()
  
  "CromwellAPIClient" should "fail the Future if the HttpResponse is unsuccessful" in {
    val errorMessage =
      """|{
         |  "status": "fail",
         |  "message": "Cromwell service shutting down"
         |}
         |""".stripMargin.trim
    val client = new CromwellClient(new URL("http://fakeurl"), "v1") {
      override def executeRequest(request: HttpRequest, headers: List[HttpHeader]): Future[HttpResponse] = Future.successful(
        new HttpResponse(
          StatusCodes.ServiceUnavailable,
          List.empty[HttpHeader],
          HttpEntity(ContentTypes.`application/json`, errorMessage),
          HttpProtocols.`HTTP/1.1`
        )
      )
    }

    for {
      exception <- recoverToExceptionIf[UnsuccessfulRequestException] {
        client.version.asIo.unsafeToFuture()
      }
      entity = exception.httpResponse.entity
      content <- Unmarshal(entity).to[String]
      _ = entity.contentType should be(ContentTypes.`application/json`)
      _ = content should be(errorMessage)
    } yield succeed
  }
}
