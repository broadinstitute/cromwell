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

  it should "emit an Http response on failure port if request is not successful and not retryable" in {
    val mockResponse = HttpResponse(status = StatusCodes.BadRequest)
    
    val httpMock = new HttpMock[Any](MockHttpResponse(Success(mockResponse), 1)).httpMock()
    val httpRetryMock = HttpFlowWithRetry[Any](httpMock)

    val (graph, probeSuccess, probeFailure) = makeRunnableGraph(httpRetryMock)

    graph.run()

    probeSuccess.expectNoMsg(1 second)
    probeFailure.expectMsgPF(1 second) {
      case (response: HttpResponse, _) => response.status shouldBe StatusCodes.BadRequest
    }

    stopProbe(probeSuccess)
    stopProbe(probeFailure)
  }

  it should "retry a request if it fails" in {
    val mockFailure = Failure[HttpResponse](new Exception("Request failed - part of test flow"))
    val mockOK = HttpResponse(status = StatusCodes.OK)
    
    val httpMock = new HttpMock[Any](MockHttpResponse(mockFailure, 1), MockHttpResponse(Success(mockOK), 1)).httpMock()
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
  
  it should "retry a request if the response is not a success and is retryable" in {
    def isRetryable(httpResponse: HttpResponse) = httpResponse.status == StatusCodes.Unauthorized
    
    val mockFailure = HttpResponse(status = StatusCodes.Unauthorized)
    val mockOK = HttpResponse(status = StatusCodes.OK)
    val httpMock = new HttpMock[Any](MockHttpResponse(Success(mockFailure), 1), MockHttpResponse(Success(mockOK), 1)).httpMock()
    val httpRetryMock = HttpFlowWithRetry[Any](httpMock, isRetryable = isRetryable)

    val (graph, probeSuccess, probeFailure) = makeRunnableGraph(httpRetryMock)

    graph.run()

    probeFailure.expectNoMsg(1 second)
    probeSuccess.expectMsgPF(1 second) {
      case (response: HttpResponse, _) => response.status shouldBe StatusCodes.OK
    }

    stopProbe(probeSuccess)
    stopProbe(probeFailure)
  }

  it should "handle flaky connections with a sufficiently large retry buffer" in {
    val mock200 = HttpResponse(status = StatusCodes.OK)
    val mock429 = HttpResponse(status = StatusCodes.TooManyRequests)
    
    // only 20% of requests come back successful, rest are out of quota errors
    val mock200s = MockHttpResponse(Success(mock200), 1)
    val mock429s = MockHttpResponse(Success(mock429), 5)

    val httpMock = new HttpMock[Any](mock200s, mock429s).httpMock()
    // Give a retry buffer of 1000 (= queues up to 1000 retry requests before starting dropping them) 
    val httpRetryMock = new HttpFlowWithRetry[Any](httpMock, retryBufferSize = 1000)

    // Send 1000 requests
    val nbMessagesToSend = 1000
    val mockRequest = HttpRequest()
    val sourceIterator = 1 to nbMessagesToSend map { _ => (mockRequest, ()) } toIterator
    val source = Source.fromIterator { () => sourceIterator }
    val (graph, probeSuccess, probeFailure) = makeRunnableGraph(httpRetryMock, source)

    graph.run()

    probeFailure.expectNoMsg(5 seconds)
    var successes: Int = 0
    
    probeSuccess.receiveWhile(5 seconds, 500 millis, nbMessagesToSend) {
      case _ => successes += 1
    }

    // All responses should come back successfully
    successes == 1000 shouldBe true

    stopProbe(probeSuccess)
    stopProbe(probeFailure)
  }

  it should "not deadlock even under high retry pressure" in {
    val mock200 = HttpResponse(status = StatusCodes.OK)
    val mock429 = HttpResponse(status = StatusCodes.TooManyRequests)

    // only 20% of requests come back successful, rest are out of quota errors
    val mock200s = MockHttpResponse(Success(mock200), 1)
    val mock429s = MockHttpResponse(Success(mock429), 5)

    val httpMock = new HttpMock[Any](mock200s, mock429s).httpMock()
    // Give a retry buffer of 10 (= queues up to 10 retry requests before starting dropping them)
    val httpRetryMock = new HttpFlowWithRetry[Any](httpMock, retryBufferSize = 10)

    // Send 1000 requests, the retry buffer will overflow, and start dropping requests
    // Note that this is an extreme case, 80% failure rate with 100 * retry buffer size elements and no throttling
    val nbMessagesToSend = 1000
    val mockRequest = HttpRequest()
    val sourceIterator = 1 to nbMessagesToSend map { _ => (mockRequest, ()) } toIterator
    val source = Source.fromIterator { () => sourceIterator }
    val (graph, probeSuccess, probeFailure) = makeRunnableGraph(httpRetryMock, source)

    graph.run()

    probeFailure.expectNoMsg(5 seconds)
    var successes: Int = 0
    // The second parameter makes sure that the stream keeps emitting responses at least every 500 milliseconds
    probeSuccess.receiveWhile(5 seconds, 500 millis, nbMessagesToSend) {
      case _ => successes += 1
    }

    // A minimum amount of requests should come back
    successes > 200 shouldBe true

    stopProbe(probeSuccess)
    stopProbe(probeFailure)
  }
}
