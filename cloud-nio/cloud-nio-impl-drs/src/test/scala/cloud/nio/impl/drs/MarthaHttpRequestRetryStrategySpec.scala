package cloud.nio.impl.drs

import java.io.IOException

import org.apache.http.HttpVersion
import org.apache.http.client.methods.CloseableHttpResponse
import org.apache.http.message.BasicStatusLine
import org.apache.http.protocol.HttpContext
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito

import scala.concurrent.duration._

class MarthaHttpRequestRetryStrategySpec extends AnyFlatSpec with Matchers with Mockito {

  behavior of "MarthaHttpRequestRetryStrategy"

  it should "retry 500 errors a configured number of times" in {
    val drsConfig = MockDrsPaths.mockDrsConfig.copy(numRetries = 3)
    val retryStrategy = new MarthaHttpRequestRetryStrategy(drsConfig)
    val http500Response = mock[CloseableHttpResponse].smart
    http500Response.getStatusLine returns new BasicStatusLine(HttpVersion.HTTP_1_1, 500, "Testing 500")
    val httpContext = mock[HttpContext].smart

    // initial failure
    retryStrategy.retryRequest(http500Response, 1, httpContext) should be(true)
    // three retries
    retryStrategy.retryRequest(http500Response, 2, httpContext) should be(true)
    retryStrategy.retryRequest(http500Response, 3, httpContext) should be(true)
    retryStrategy.retryRequest(http500Response, 4, httpContext) should be(true)
    // no more retries
    retryStrategy.retryRequest(http500Response, 5, httpContext) should be(false)
  }

  it should "retry 500 errors even after a number of 429 errors" in {
    val drsConfig = MockDrsPaths.mockDrsConfig.copy(numRetries = 3)
    val retryStrategy = new MarthaHttpRequestRetryStrategy(drsConfig)
    val http500Response = mock[CloseableHttpResponse].smart
    http500Response.getStatusLine returns new BasicStatusLine(HttpVersion.HTTP_1_1, 500, "Testing 500")
    val http429Response = mock[CloseableHttpResponse].smart
    http429Response.getStatusLine returns new BasicStatusLine(HttpVersion.HTTP_1_1, 429, "Testing 429")
    val httpContext = mock[HttpContext].smart

    // initial failure
    retryStrategy.retryRequest(http500Response, 1, httpContext) should be(true)
    // one retry
    retryStrategy.retryRequest(http500Response, 2, httpContext) should be(true)
    // a couple 429s
    retryStrategy.retryRequest(http429Response, 3, httpContext) should be(true)
    retryStrategy.retryRequest(http429Response, 4, httpContext) should be(true)
    // two more 500s should still retry
    retryStrategy.retryRequest(http500Response, 5, httpContext) should be(true)
    retryStrategy.retryRequest(http500Response, 6, httpContext) should be(true)
    // can still retry 429
    retryStrategy.retryRequest(http429Response, 7, httpContext) should be(true)
    retryStrategy.retryRequest(http429Response, 8, httpContext) should be(true)
    // but no more retries for 500
    retryStrategy.retryRequest(http500Response, 9, httpContext) should be(false)
  }

  it should "not retry an HTTP 401" in {
    val drsConfig = MockDrsPaths.mockDrsConfig.copy(numRetries = 3)
    val retryStrategy = new MarthaHttpRequestRetryStrategy(drsConfig)
    val http400Response = mock[CloseableHttpResponse].smart
    http400Response.getStatusLine returns new BasicStatusLine(HttpVersion.HTTP_1_1, 401, "Testing 401")
    val httpContext = mock[HttpContext].smart

    retryStrategy.retryRequest(http400Response, 1, httpContext) should be(false)
  }

  it should "retry IO exceptions a configured number of times" in {
    val drsConfig = MockDrsPaths.mockDrsConfig.copy(numRetries = 3)
    val retryStrategy = new MarthaHttpRequestRetryStrategy(drsConfig)
    val exception = mock[IOException].smart
    val httpContext = mock[HttpContext].smart

    // initial failure
    retryStrategy.retryRequest(exception, 1, httpContext) should be(true)
    // three retries
    retryStrategy.retryRequest(exception, 2, httpContext) should be(true)
    retryStrategy.retryRequest(exception, 3, httpContext) should be(true)
    retryStrategy.retryRequest(exception, 4, httpContext) should be(true)
    // no more retries
    retryStrategy.retryRequest(exception, 5, httpContext) should be(false)
  }

  it should "retry only up to a configured maximum duration" in {
    val drsConfig =
      MockDrsPaths.mockDrsConfig.copy(
        waitInitial = 10.seconds,
        waitMultiplier = 2.0d,
        waitMaximum = 1.minute,
        waitRandomizationFactor = 0d,
      )
    val retryStrategy = new MarthaHttpRequestRetryStrategy(drsConfig)

    retryStrategy.getRetryInterval should be(10000L)
    retryStrategy.getRetryInterval should be(20000L)
    retryStrategy.getRetryInterval should be(40000L)
    retryStrategy.getRetryInterval should be(60000L)
    retryStrategy.getRetryInterval should be(60000L)
  }
}
