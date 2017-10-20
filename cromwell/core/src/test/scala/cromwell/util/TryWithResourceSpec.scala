package cromwell.util

import org.scalatest.{FlatSpec, Matchers}
import TryWithResource._

import scala.util.{Failure, Success}

class TryWithResourceSpec extends FlatSpec with Matchers {
  behavior of "tryWithResource"

  it should "catch instantiation errors" in {
    val triedMyBest = tryWithResource(() => if (1 == 1) throw InstantiationException else null) { _ => 5 }
    triedMyBest should be(Failure(InstantiationException))
  }

  it should "close the closeable" in {
    val myCloseable = new MyCloseable
    val triedMyBest = tryWithResource(() => myCloseable) { _.value } // Nothing special about 5... Just need to return something!
    triedMyBest should be(Success(5))
    myCloseable.isClosed should be(true)
  }

  it should "catch errors and still close the closeable" in {
    val myCloseable = new MyCloseable
    val triedMyBest = tryWithResource(() => myCloseable) { _.badValue }
    triedMyBest should be(Failure(ReadValueException))
    myCloseable.isClosed should be(true)
  }

  it should "be robust to failures in close methods" in {
    val myCloseable = new FailingCloseable
    val triedMyBest = tryWithResource(() => myCloseable) { _.value }
    triedMyBest should be(Failure(CloseCloseableException))
    val triedMyBest2 = tryWithResource(() => myCloseable) { _.badValue }
    triedMyBest2 match {
      case Failure(ReadValueException) => ReadValueException.getSuppressed.headOption should be(Some(CloseCloseableException))
      case x => fail(s"$x was not equal to $ReadValueException")
    }
  }
}

class MyCloseable extends AutoCloseable {
  var isClosed = false

  val value = if (isClosed) throw ReadValueException else 5 // Ensures we aren't closed when .value is called
  def badValue = throw ReadValueException

  override def close() = {
    isClosed = true
  }
}

class FailingCloseable extends MyCloseable {
  override def close(): Unit = throw CloseCloseableException
}

case object InstantiationException extends Exception("Oh Teh Noes instantiating!")
case object ReadValueException extends Exception("Oh Teh Noes reading!")
case object CloseCloseableException extends Exception("Oh Teh Noes closing!")
