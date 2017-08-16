package cromwell.docker.registryv2.flows

import akka.NotUsed
import akka.actor.Scheduler
import akka.http.scaladsl.model.{HttpRequest, HttpResponse, StatusCodes}
import akka.stream.{ActorMaterializer, FanOutShape2}
import akka.stream.javadsl.MergePreferred
import akka.stream.scaladsl.{Flow, GraphDSL, Partition}
import cromwell.docker.registryv2.flows.FlowUtils._
import cromwell.docker.registryv2.flows.HttpFlowWithRetry._
import cromwell.core.retry.{Backoff, SimpleExponentialBackoff}
import org.slf4j.LoggerFactory

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object HttpFlowWithRetry {
  val Logger = LoggerFactory.getLogger("HttpLogger")
  
  def isRetryable(response: HttpResponse) = {
    response.status match {
      case StatusCodes.InternalServerError => true
      case StatusCodes.ServiceUnavailable => true
      case StatusCodes.BadGateway => true
      case StatusCodes.GatewayTimeout => true
      case StatusCodes.RequestTimeout => true
      case other @ _ => isTransient(response)
    }
  }

  def isTransient(response: HttpResponse) = {
    response.status match {
      case StatusCodes.TooManyRequests => true
      case _ => false
    }
  }

  /**
    * In order to allow for retries, the http context object needs to encapsulate the original request,
    * so that it can be re-submitted if necessary.
    * This type provides this encapsulation by creating a pair of (T, HttpRequest) where T is the custom type that
    * the user wishes to use as a context for their specific use case.
    */
  case class ContextWithRequest[T](userContext: T, request: HttpRequest, currentAttempt: Int, backoff: Backoff) {
    // Back off is mutable ! This will change every time !
    def retryIn = backoff.backoffMillis.millis
    
    def withNextAttempt = copy(currentAttempt = currentAttempt + 1)
  }

  /**
    * Type for a http flow using ContextWithRequest as context type.
    */
  type RetryableHttpFlow[T] = Flow[(HttpRequest, ContextWithRequest[T]), (Try[HttpResponse], ContextWithRequest[T]), NotUsed]
}

/**
  * Adds additional retry logic around an existing akka http flow.
  * @param httpClientFlow The akka http flow to use underneath
  * @param retryBufferSize size of the buffer for requests to be retried. When the buffer is full, requests WILL BE DROPPED
  *                        Act accordingly by resubmitting requests after a certain amount of time with no response for example.
  * @tparam T Type of the context data to be passed along with the request
  * @return A http flow that will retry on request failures and on failed responses as seen fit by the isRetryable method
  */
case class HttpFlowWithRetry[T](
                            httpClientFlow: RetryableHttpFlow[T],
                            retryBufferSize: Int = 100,
                            requestBackoff: () => Backoff = () => SimpleExponentialBackoff(1 second, 2 minutes, 3D),
                            maxAttempts: Int = 3
                          )(implicit val scheduler: Scheduler, ec: ExecutionContext, mat: ActorMaterializer) {
  
  lazy val flow = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    // Http shape
    val http = builder.add(httpClientFlow)
    
    // Wrap the user context into a pair of (flowContext, httpRequest) to allow for retries
    val source = builder.add(Flow[(HttpRequest, T)] map {
      case (httpRequest, flowContext) => (httpRequest, ContextWithRequest(flowContext, httpRequest, 1, requestBackoff()))
    })
    
    // Partition Try[HttpResponse] into Success and Failure
    val partitionResponseTry = builder.add(fanOutTry[HttpResponse, ContextWithRequest[T]])
    
    // Those are pairs like: (HttpResponse, Context) 
    val requestSuccessful = partitionResponseTry.out0
    
    // Those are pairs like: (Throwable, Context) 
    // -> akka-http already retries request if they fail, so we won't retry this and send it to the failure port
    val requestFailed = partitionResponseTry.out1 map {
      case (failure, ContextWithRequest(flowContext, _, _, _)) => failure -> flowContext
    }

    // Splits successful http response into 2 categories:
    // - out(0) emits terminal http resonse (successful or unsuccessful non-retryable)
    // - out(1) emits unsusccessful but retryable responses
    val partitionHttpResponse = builder.add(Partition[(HttpResponse, ContextWithRequest[T])](2, {
      // Successful return code
      case (httpResponse, flowContext) if httpResponse.status.isSuccess() || !shouldRetry(httpResponse, flowContext) => 0
      // Failed return code but retryable
      case (_, _) => 1
    }))
    
    // Merges requests coming from 3 input ports:
    // - in(0) takes requests from the outside
    // - in(1) takes retryable failed responses (with non successful and retryable status code)
    val mergeRequests = builder.add(MergePreferred.create[(HttpRequest, ContextWithRequest[T])](1))

    // Either success or non retryable unsuccessful http response, either way it won't be retried
    val terminalResponse = partitionHttpResponse.out(0) map {
      case (response, ContextWithRequest(flowContext, _, _, _)) => response -> flowContext
    }
    
    // Retryable http response
    val retryableResponse = partitionHttpResponse.out(1)
      // Create a new request, yields a future that will complete after the appropriate backoff time
      .map(toRetryableRequest)
      // wait for the future to complete, retryBufferSize sets how many retryable future we can wait for in parallel
      .mapAsync(retryBufferSize)(identity)
    
    // Outside requests coming in
    source ~> mergeRequests.in(0)
    
    // Submit request to underlying http flow  -  Partition responses: Failure | Success
    mergeRequests     ~>     http      ~>         partitionResponseTry.in
    
                                                  // Success -> analyze response and partition into 2 (see above)
                                                  requestSuccessful ~> partitionHttpResponse
    
    // Retryable http response are sent back for retry
    mergeRequests.preferred          <~           retryableResponse.outlet

    new FanOutShape2(source.in, terminalResponse.outlet, requestFailed.outlet)
  }
  
  private def shouldRetry(httpResponse: HttpResponse, contextWithRequest2: ContextWithRequest[T]) = {
    isRetryable(httpResponse) && contextWithRequest2.currentAttempt <= maxAttempts
  }

  /**
    * Create a re-submittable request from a failed retryable
    * @return a future that will complete after the appropriate backoff time
    */
  private def toRetryableRequest(value: (HttpResponse, ContextWithRequest[T])) = value match {
    case (response, contextWithRequest) => 
      val nextRetryIn = contextWithRequest.retryIn
      val nextRequest = (contextWithRequest.request, contextWithRequest.withNextAttempt)
      // This response will never be consumed by anyone, so discard its content here to avoid pool freeze
      // http://doc.akka.io/docs/akka-http/10.0.5/scala/http/client-side/request-level.html#using-the-future-based-api-in-actors
      // https://github.com/akka/akka/issues/19538
      // https://github.com/akka/akka-http/issues/183
      // https://github.com/akka/akka-http/issues/117
      response.discardEntityBytes().future() flatMap { _ =>
        akka.pattern.after(nextRetryIn, scheduler) {
          Future.successful(nextRequest)
        }
      } recover {
        case failure =>
          // Can't do much here except log the error and keep going with the next request
          Logger.error(s"Failed to discard entity bytes for response $response", failure)
          nextRequest
      }
  }
}
