package cromwell.backend.impl.aws.util

import com.amazonaws.AmazonWebServiceRequest
import com.amazonaws.handlers.AsyncHandler
import cromwell.backend.impl.aws.util.AwsSdkAsyncHandler._

import scala.concurrent.Promise

/**
  * Scala wrapper around the AWS SDK's AsyncHandler callbacks
  */
final class AwsSdkAsyncHandler[REQUEST <: AmazonWebServiceRequest, RESULT] extends AsyncHandler[REQUEST, RESULT] {
  private val promise = Promise[AwsSdkAsyncResult[REQUEST, RESULT]]()
  override def onError(exception: Exception): Unit = {
    promise.tryFailure(exception)
    ()
  }
  override def onSuccess(request: REQUEST, result: RESULT): Unit = {
    promise.trySuccess(AwsSdkAsyncResult(request, result))
    ()
  }

  def future = promise.future
}

object AwsSdkAsyncHandler {
  case class AwsSdkAsyncResult[REQUEST <: AmazonWebServiceRequest, RESULT](request: REQUEST, result: RESULT)
}
