package cromwell.core.callcaching.docker.registryv2.flows

import akka.NotUsed
import akka.http.scaladsl.model.{HttpRequest, HttpResponse, StatusCodes}
import akka.stream.scaladsl.{Flow, GraphDSL, Merge, Partition}
import akka.stream.{FanOutShape2, OverflowStrategy}
import cromwell.core.callcaching.docker.registryv2.flows.FlowUtils._
import cromwell.core.callcaching.docker.registryv2.flows.HttpFlowWithRetry._

import scala.util.Try

object HttpFlowWithRetry {
  def defaultIsRetryable(response: HttpResponse) = {
    response.status match {
      case StatusCodes.TooManyRequests => true
      case StatusCodes.InternalServerError => true
      case _ => false
    }
  }

  /**
    * In order to allow for retries, the http context object needs to encapsulate the original request,
    * so that it can be re-submitted if necessary.
    * This type provides this encapsulation by creating a pair of (T, HttpRequest) where T is the custom type that
    * the user wishes to use as a context for his specific use case.
    */
  type ContextWithRequest[T] = (T, HttpRequest)

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
  * @param isRetryable function to determine weather or not an http request should be retried based on the http response                        
  * @tparam T Type of the context data to be passed along with the request
  * @return A http flow that will retry on request failures and on failed responses as seen fit by the isRetryable method
  */
case class HttpFlowWithRetry[T](
                            httpClientFlow: RetryableHttpFlow[T],
                            retryBufferSize: Int = 1000,
                            isRetryable: HttpResponse => Boolean = defaultIsRetryable
                          ) {

  lazy val flow = GraphDSL.create() { implicit builder =>
    import GraphDSL.Implicits._

    // Http shape
    val http = builder.add(httpClientFlow)
    
    // Wrap the user context into a pair of (flowContext, httpRequest) to allow for retries
    val source = builder.add(Flow[(HttpRequest, T)] map {
      case (httpRequest, flowContext) => (httpRequest, (flowContext, httpRequest))
    })
    
    // Partition Try[HttpResponse] into Success and Failure
    val partitionResponseTry = builder.add(fanOutTry[HttpResponse, (T, HttpRequest)])
    
    val requestSuccessful = partitionResponseTry.out0
    // Extract the request from the failed try so it can be resubmitted
    val requestFailed = partitionResponseTry.out1.log("request failed") map toRetryableRequest
    
    // Merges requests coming from 3 input ports:
    // - in(0) takes requests from the outside
    // - in(1) takes retryable failed requests
    // - in(2) takes retryable failed responses (with non successful and retryable status code)
    val mergeRequests = builder.add(Merge[(HttpRequest, (T, HttpRequest))](3))
    
    // Splits http response into 3 categories:
    // - out(0) emits successful responses
    // - out(1) emits failed but retryable responses
    // - out(2) emits failed, non-retryable responses
    val partitionHttpResponse = builder.add(Partition[(HttpResponse, (T, HttpRequest))](3, {
      // Successful return code
      case (httpResponse, flowContext) if httpResponse.status.isSuccess() => 0
      // Failed return code but retryable
      case (httpResponse, flowContext) if isRetryable(httpResponse) => 1
      // Failed return code an non retryable
      case (httpResponse, flowContext) => 2
    }))

    val responseSuccessful = partitionHttpResponse.out(0) map toTerminalResponse
    val responseRetryable = partitionHttpResponse.out(1) map toRetryableRequest
    val responseFailed = partitionHttpResponse.out(2) map toTerminalResponse
    
    // Buffer for retryable requests. Will DROP requests when full
    // Will also delay re submission
    val retryBuffer = Flow[(HttpRequest, (T, HttpRequest))]
      .buffer(retryBufferSize, OverflowStrategy.dropHead)

    // Submit outside requests to the first port of the mergeRequests merger 
    source ~> mergeRequests.in(0)
    
    // Submit request to underlying http flow  -  Partition responses: Failure | Success
    mergeRequests     ~>     http      ~>         partitionResponseTry.in
                                                        // Success -> analyze response and partition into 3 (see above)
                                                        requestSuccessful ~> partitionHttpResponse
    // Retry failed requests (Try[HttpResponse] was a Failure)
    mergeRequests.in(1)     <~ retryBuffer <~     requestFailed.outlet
    // Retry retryable failed responses (the HttpResponse had a non successful return code and was deeemed retryable)
    mergeRequests.in(2)        <~        retryBuffer         <~        responseRetryable.outlet

    new FanOutShape2(source.in, responseSuccessful.outlet, responseFailed.outlet)
  }

  // create a re-submittable request from a failed retryable
  private def toRetryableRequest(value: (Any, (T, HttpRequest))) = value match {
    case (_, (flowContext, request)) => (request, (flowContext, request))
  }

  // extract only the response and the context so it can be emitted through the output port
  private def toTerminalResponse(value: (HttpResponse, (T, HttpRequest))) = value match {
    case (response, (flowContext, _)) => (response, flowContext)
  }
}
