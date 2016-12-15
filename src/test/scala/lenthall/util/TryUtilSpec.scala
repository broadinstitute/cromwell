package lenthall.util

import lenthall.exception.AggregatedException
import lenthall.util.TryUtil._
import org.scalatest.{FlatSpec, Matchers}

import scala.compat.Platform.EOL
import scala.util.{Failure, Success}

class TryUtilSpec extends FlatSpec with Matchers {

  behavior of "TryUtil"

  it should "not stringify successes" in {
    val result = stringifyFailures(Traversable(Success("success")))
    result should be(empty)
  }

  it should "stringify failures" in {
    val result = stringifyFailures(Traversable(Failure(new RuntimeException("failed"))))
    result should have size 1
    result.head should startWith(s"java.lang.RuntimeException: failed$EOL")
  }

  it should "sequence successful seqs" in {
    val result = sequence(Seq(Success("success")), "prefix")
    result.isSuccess should be(true)
    result.get should contain theSameElementsAs Seq("success")
  }

  it should "sequence failed seqs" in {
    val result = sequence(Seq(Failure(new RuntimeException("failed"))), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 1
    exception.throwables.head.getMessage should be("failed")
  }

  it should "sequence successful maps" in {
    val result = sequenceMap(Map("key" -> Success("success")), "prefix")
    result.isSuccess should be(true)
    result.get should contain theSameElementsAs Map("key" -> "success")
  }

  it should "sequence failed maps" in {
    val result = sequenceMap(Map("key" -> Failure(new RuntimeException("failed"))), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 1
    exception.throwables.head.getMessage should be("failed")
  }

}
