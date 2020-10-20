package cloud.nio.impl.drs

import java.io.IOException

import org.apache.http.HttpResponse
import org.apache.http.client.{HttpRequestRetryHandler, ServiceUnavailableRetryStrategy}
import org.apache.http.protocol.HttpContext

import scala.concurrent.duration._

class MarthaHttpRequestRetryStrategy(drsConfig: DrsConfig)
  extends ServiceUnavailableRetryStrategy with HttpRequestRetryHandler {

  // We can execute a total of one time, plus the number of retries
  private val executionMax: Int = drsConfig.numRetries + 1
  private val waitMultiplier: Double = drsConfig.waitMultiplier
  private val waitMax: FiniteDuration = drsConfig.waitMaximum
  private var waitNext: Duration = drsConfig.waitInitial min waitMax
  private var transientFailures: Int = 0

  /** Returns true if an IOException should be immediately retried. */
  override def retryRequest(exception: IOException, executionCount: Int, context: HttpContext): Boolean = {
    retryRequest(executionCount)
  }

  /** Returns true if HttpResponse should be retried after getRetryInterval. */
  override def retryRequest(response: HttpResponse, executionCount: Int, context: HttpContext): Boolean = {
    response.getStatusLine.getStatusCode match {
      case 429 => retryRequestTransient(executionCount)
      case code if 500 <= code && code <= 599 => retryRequest(executionCount)
      case _ => false
    }
  }

  /** Returns the number of milliseconds to wait before retrying an HttpResponse. */
  override def getRetryInterval: Long = {
    val waitCurrent = waitNext
    waitNext = (waitNext * waitMultiplier) min waitMax
    waitCurrent.toMillis
  }

  private def retryRequestTransient(executionCount: Int): Boolean = {
    transientFailures += 1
    retryRequest(executionCount)
  }

  private def retryRequest(executionCount: Int): Boolean = {
    // The first execution is executionCount == 1
    executionCount - transientFailures <= executionMax
  }
}
