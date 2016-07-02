package cromwell.util

import cromwell.core.CromwellFatalException
import cromwell.core.retry.SimpleExponentialBackoff
import org.scalatest.mock.MockitoSugar
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.Success

class TryUtilSpec extends FlatSpec with Matchers with MockitoSugar {

  class TransientException extends Exception
  class MockWork {
    var counter: Int = _

    /**
      * @param n number of times an exception is thrown before succeeding
      * @param transients how many transient exceptions to raise (must be <= n)
      */
    def failNTimes(n: Int, transients: Int = 0): Option[Int] => Int = {
      counter = n
      def func(prior: Option[Int]): Int = {
        if (counter > 0) {
          counter -= 1
          if (counter <= transients) throw new TransientException
          else throw new IllegalArgumentException("Failed")
        }
        9
      }
      func
    }
  }

  //val logger = mock[WorkflowLogger]
  val backoff = SimpleExponentialBackoff(50 milliseconds, 10 seconds, 1D)
}
