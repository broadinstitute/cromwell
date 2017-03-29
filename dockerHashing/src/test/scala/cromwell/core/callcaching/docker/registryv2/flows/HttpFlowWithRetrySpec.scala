package cromwell.core.callcaching.docker.registryv2.flows

import akka.NotUsed
import akka.http.scaladsl.model.{HttpRequest, HttpResponse, StatusCodes}
import akka.stream._
import akka.stream.scaladsl.{GraphDSL, RunnableGraph, Sink, Source}
import akka.testkit.{ImplicitSender, TestProbe}
import cromwell.core.TestKitSuite
import cromwell.core.callcaching.docker.{HttpMock, MockHttpResponse}
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}

class HttpFlowWithRetrySpec extends TestKitSuite with FlatSpecLike with Matchers with ImplicitSender {
  behavior of "HttpFlowWithRetry"

  implicit val materializer = ActorMaterializer()
  implicit val scheduler = system.scheduler
  implicit val ec = system.dispatcher

  override protected def afterAll() = {
    materializer.shutdown()
    super.afterAll()
  }
  
  def stopProbe(probe: TestProbe) = system stop probe.ref

  def makeRunnableGraph(httpFlowWithRetry: HttpFlowWithRetry[Any],
                        source: Source[(HttpRequest, Unit), NotUsed] = Source.single((HttpRequest(), ()))
                       ) = {
    
    val probeSuccess = TestProbe()
    val probeFailure = TestProbe()
    val sinkSuccess = Sink.actorRef(probeSuccess.ref, null)
    val sinkFailure = Sink.actorRef(probeFailure.ref, null)
    
    val graph = RunnableGraph.fromGraph(GraphDSL.create() { implicit b =>
      import akka.stream.scaladsl.GraphDSL.Implicits._
      val httpRetryMockPorts = b.add(httpFlowWithRetry.flow)

      source ~> httpRetryMockPorts.in
      httpRetryMockPorts.out0 ~> sinkSuccess
      httpRetryMockPorts.out1 ~> sinkFailure

      ClosedShape
    })
    
    (graph, probeSuccess, probeFailure)
  }
  
  it should "emit an Http response on success port if request is successful" in {
    val mockResponse = HttpResponse(status = StatusCodes.OK)
    
    val httpMock = new HttpMock[Any](MockHttpResponse(Success(mockResponse), 1)).httpMock()
    val httpRetryMock = HttpFlowWithRetry[Any](httpMock)
    
    val (graph, probeSuccess, probeFailure) = makeRunnableGraph(httpRetryMock)

    graph.run()

    probeFailure.expectNoMsg(1 second)
    probeSuccess.expectMsgPF(1 second) {
      case (response: HttpResponse, _) => response.status shouldBe StatusCodes.OK
    }

    stopProbe(probeSuccess)
    stopProbe(probeFailure)
  }

  it should "send back the http response if it is not successful and not retryable" in {
    val mockResponse = HttpResponse(status = StatusCodes.BadRequest)
    
    val httpMock = new HttpMock[Any](MockHttpResponse(Success(mockResponse), 1)).httpMock()
    val httpRetryMock = HttpFlowWithRetry[Any](httpMock)

    val (graph, probeSuccess, probeFailure) = makeRunnableGraph(httpRetryMock)

    graph.run()

    probeFailure.expectNoMsg(1 second)
    probeSuccess.expectMsgPF(1 second) {
      case (response: HttpResponse, _) => response.status shouldBe StatusCodes.BadRequest
    }

    stopProbe(probeSuccess)
    stopProbe(probeFailure)
  }

  it should "fail a request if it fails" in {
    val mockFailure = Failure[HttpResponse](new Exception("Request failed - part of test flow"))
    val mockOK = HttpResponse(status = StatusCodes.OK)
    
    val httpMock = new HttpMock[Any](MockHttpResponse(mockFailure, 1), MockHttpResponse(Success(mockOK), 1)).httpMock()
    val httpRetryMock = HttpFlowWithRetry[Any](httpMock)

    val (graph, probeSuccess, probeFailure) = makeRunnableGraph(httpRetryMock)

    graph.run()

    probeSuccess.expectNoMsg(1 second)
    probeFailure.expectMsgPF(1 second) {
      case (failure: Throwable, _) => failure shouldBe mockFailure.failed.get
    }

    stopProbe(probeSuccess)
    stopProbe(probeFailure)
  }
  
  it should "retry a request if the response is not a success and is retryable" in {
    val mockFailure = HttpResponse(status = StatusCodes.RequestTimeout)
    val mockOK = HttpResponse(status = StatusCodes.OK)
    val httpMock = new HttpMock[Any](MockHttpResponse(Success(mockFailure), 1), MockHttpResponse(Success(mockOK), 1)).httpMock()
    val httpRetryMock = HttpFlowWithRetry[Any](httpMock)

    val (graph, probeSuccess, probeFailure) = makeRunnableGraph(httpRetryMock)

    graph.run()

    probeFailure.expectNoMsg(1 second)
    probeSuccess.expectMsgPF(1 second) {
      case (response: HttpResponse, _) => response.status shouldBe StatusCodes.OK
    }

    stopProbe(probeSuccess)
    stopProbe(probeFailure)
  }

  it should "exponentially backoff and succeed if response is successful before max attempts" in {
    val mock200 = HttpResponse(status = StatusCodes.OK)
    val mock408 = HttpResponse(status = StatusCodes.RequestTimeout)

    val mock408s = MockHttpResponse(Success(mock408), 3)
    val mock200s = MockHttpResponse(Success(mock200), 1)

    val httpMock = new HttpMock[Any](mock408s, mock200s).httpMock()
    // Give a retry buffer of 10 (= queues up to 10 retry requests before starting dropping them)
    val httpRetryMock = new HttpFlowWithRetry[Any](httpMock, retryBufferSize = 10)
    
    val (graph, probeSuccess, probeFailure) = makeRunnableGraph(httpRetryMock)

    graph.run()

    probeFailure.expectNoMsg(1 second)
    // backoff starts at 1 second, multiplier is 3 and default randomization is 20%
    // which means, worst case scenario, we need to wait
    // (1 + 0.2) + (3 + 0.6) + (9 + 1.8) = 15.06 -> 20 to account for test latency
    probeSuccess.expectMsgPF(20 seconds) {
      case (response: HttpResponse, _) => response.status shouldBe StatusCodes.OK
    }

    stopProbe(probeSuccess)
    stopProbe(probeFailure)
  }

  it should "exponentially backoff and fail if response is still unsuccessful after max attempts" in {
    val mock408 = HttpResponse(status = StatusCodes.RequestTimeout)

    val mock408s = MockHttpResponse(Success(mock408), 1)

    val httpMock = new HttpMock[Any](mock408s).httpMock()
    // Give a retry buffer of 10 (= queues up to 10 retry requests before starting dropping them)
    val httpRetryMock = new HttpFlowWithRetry[Any](httpMock, retryBufferSize = 10)

    val (graph, probeSuccess, probeFailure) = makeRunnableGraph(httpRetryMock)

    graph.run()

    probeFailure.expectNoMsg(1 second)
    probeSuccess.expectMsgPF(20 seconds) {
      case (response: HttpResponse, _) => response.status shouldBe StatusCodes.RequestTimeout
    }

    stopProbe(probeSuccess)
    stopProbe(probeFailure)
  }
}
