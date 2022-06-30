package cromwell.filesystems.oss.nio

import cromwell.core.TestKitSuite

import scala.util.{Success, Failure, Try}

case class FatalError(message: String = "", cause: Throwable = None.orNull) extends Exception(message, cause)
case class TransientError(message: String = "", cause: Throwable = None.orNull) extends Exception(message, cause)
case class RetryableError(message: String = "", cause: Throwable = None.orNull) extends Exception(message, cause)

class RetryContext {
  var retried = 0

  def doSth(f: Int => Try[Int]): Unit = {
    f(retried) match {
      case Success(_) => retried += 1
      case Failure(e: RetryableError) =>
        retried += 1
        throw e
      case Failure(e: TransientError) =>
        retried += 1
        throw e
      case Failure(e) => throw e
    }
  }
}

object OssStorageRetrySpec {
  def isFatal(t: Throwable): Boolean = t match {
    case _: FatalError => true
    case _ => false
  }

  def isTransient(t: Throwable): Boolean = t match {
    case _: TransientError => true
    case _ => false
  }
}

class OssStorageRetrySpec extends TestKitSuite with OssNioUtilSpec {

  import OssStorageRetrySpec._

  behavior of s"OssStorageRetrySpec"

  it should "retry throw immediately when fatal error occours" in {
    val f = (x: Int) => if (x == 0) Failure(new FatalError) else Success(x)
    val ctx = new RetryContext()
    an [FatalError] should be thrownBy OssStorageRetry.fromTry(
      () => Try{
        ctx.doSth(f)
      },
      isFatal = isFatal,
      isTransient = isTransient
    )

    ctx.retried shouldEqual(0)
  }

  it should "retry if non-fatal error occurs" in {
    val needRetry = 5
    val f = (x: Int) => {
      if (x < needRetry) {
        Failure(new RetryableError())
      } else {
        Failure(new FatalError())
      }
    }

    val ctx = new RetryContext()
    an [FatalError] should be thrownBy OssStorageRetry.fromTry(
      () => Try{
        ctx.doSth(f)
      },
      isFatal = isFatal,
      isTransient = isTransient
    )

    ctx.retried shouldEqual(needRetry)
  }

  it should "success after retry max retries " in {

    val needRetry = 5
    val f = (x: Int) => {
      if (x < needRetry) {
        Failure(new RetryableError())
      } else {
        Success(x)
      }
    }

    val ctx = new RetryContext()
    OssStorageRetry.fromTry(
      () => Try{
        ctx.doSth(f)
      },
      isFatal = isFatal,
      isTransient = isTransient
    )

    ctx.retried shouldEqual(needRetry + 1)
  }

  it should "retry at most max retry times" in {
    val needRetry = OssStorageRetry.DEFAULT_MAX_RETRIES + 1
    val f = (x: Int) => {
      if (x < needRetry) {
        Failure(new RetryableError())
      } else {
        Success(x)
      }
    }

    val ctx = new RetryContext()
    an [RetryableError] should be thrownBy OssStorageRetry.fromTry(
      () => Try{
        ctx.doSth(f)
      },
      isFatal = isFatal,
      isTransient = isTransient
    )

    ctx.retried shouldEqual(OssStorageRetry.DEFAULT_MAX_RETRIES)
  }
}
