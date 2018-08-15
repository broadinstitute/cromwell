package cloud.nio.spi

import com.typesafe.config.Config

import scala.util.{Failure, Success, Try}
import scala.concurrent.duration._
import net.ceedubs.ficus.Ficus._

class CloudNioRetry(config: Config) {

  def defaultMaxRetries: Int = config.getAs[Int]("default-max-retries").getOrElse(0)

  def isFatal(exception: Exception): Boolean = !isTransient(exception)

  def isTransient(exception: Exception): Boolean = false

  def defaultBackOff: CloudNioBackoff = CloudNioSimpleExponentialBackoff(1.second, 60.seconds, 1.1)

  def fromTry[A](
    f: () => Try[A],
    maxRetries: Option[Int] = Option(defaultMaxRetries),
    backoff: CloudNioBackoff = defaultBackOff
  ): A = {
    val delay = backoff.backoffMillis

    f() match {
      case Success(ret)                                        => ret
      case Failure(exception: Exception) if isFatal(exception) => throw exception
      case Failure(exception: Exception) if !isFatal(exception) =>
        val retriesLeft = if (isTransient(exception)) maxRetries else maxRetries map { _ - 1 }
        if (retriesLeft.forall(_ > 0)) {
          Thread.sleep(delay)
          fromTry(f, retriesLeft, backoff.next)
        } else {
          throw exception
        }
      case Failure(throwable: Throwable) => throw throwable
    }
  }

  def from[A](f: () => A, maxRetries: Option[Int] = Option(defaultMaxRetries), backoff: CloudNioBackoff = defaultBackOff): A = {
    fromTry[A](
      () => Try(f()),
      maxRetries,
      backoff
    )
  }
}
