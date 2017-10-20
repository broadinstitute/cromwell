package lenthall.util

import lenthall.exception.AggregatedException
import lenthall.util.TryUtil._
import org.scalatest.{FlatSpec, Matchers}

import scala.compat.Platform.EOL
import scala.util.{Failure, Success, Try}

class TryUtilSpec extends FlatSpec with Matchers {

  behavior of "TryUtil"

  it should "not stringify successes" in {
    val result: Traversable[String] = stringifyFailures(Traversable(Success("success")))
    result should be(empty)
  }

  it should "stringify failures" in {
    val result: Traversable[String] = stringifyFailures(Traversable(Failure(new RuntimeException("failed"))))
    result should have size 1
    result.head should startWith(s"java.lang.RuntimeException: failed$EOL")
  }

  it should "sequence successful seqs" in {
    val result: Try[Seq[String]] = sequence(Seq(Success("success")), "prefix")
    result.isSuccess should be(true)
    result.get should contain theSameElementsAs Seq("success")
  }

  it should "sequence failed seqs" in {
    val result: Try[Seq[String]] = sequence(Seq(Failure(new RuntimeException("failed"))), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 1
    exception.throwables.head.getMessage should be("failed")
  }

  it should "sequence successful nones" in {
    val result: Try[Option[String]] = sequenceOption(None, "prefix")
    result.isSuccess should be(true)
    result.get should be(None)
  }

  it should "sequence successful somes" in {
    val result: Try[Option[String]] = sequenceOption(Option(Success("success")), "prefix")
    result.isSuccess should be(true)
    result.get should be(Option("success"))
  }

  it should "sequence failed somes" in {
    val result: Try[Option[String]] = sequenceOption(Option(Failure(new RuntimeException("failed"))), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 1
    exception.throwables.head.getMessage should be("failed")
  }

  it should "sequence successful maps" in {
    val result: Try[Map[String, String]] = sequenceMap(Map("key" -> Success("success")), "prefix")
    result.isSuccess should be(true)
    result.get should contain theSameElementsAs Map("key" -> "success")
  }

  it should "sequence failed maps" in {
    val result: Try[Map[String, String]] = sequenceMap(Map("key" -> Failure(new RuntimeException("failed"))), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 1
    exception.throwables.head.getMessage should be("failed")
  }

  it should "sequence successful keys and successful values" in {
    val result: Try[Map[String, String]] = sequenceKeyValues(
      Map(Success("success key") -> Success("success value")), "prefix")
    result.isSuccess should be(true)
    result.get should contain theSameElementsAs Map("success key" -> "success value")
  }

  it should "sequence successful keys and failed values" in {
    val result: Try[Map[String, String]] = sequenceKeyValues(
      Map(Success("success key") -> Failure(new RuntimeException("failed value"))), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 1
    exception.throwables.head.getMessage should be("failed value")
  }

  it should "sequence failed keys and successful values" in {
    val result: Try[Map[String, String]] = sequenceKeyValues(
      Map(Failure(new RuntimeException("failed key")) -> Success("success value")), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 1
    exception.throwables.head.getMessage should be("failed key")
  }

  it should "sequence failed keys and failed values" in {
    val result: Try[Map[String, String]] = sequenceKeyValues(
      Map(Failure(new RuntimeException("failed key")) -> Failure(new RuntimeException("failed value"))), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 2
    exception.throwables.head.getMessage should be("failed key")
    exception.throwables.tail.head.getMessage should be("failed value")
  }

  it should "sequence a successful key with a successful value" in {
    val result: Try[(String, String)] = sequenceTuple((Success("success key"), Success("success value")), "prefix")
    result.isSuccess should be(true)
    result.get should be(("success key", "success value"))
  }

  it should "sequence a successful key with a failed value" in {
    val result: Try[(String, String)] = sequenceTuple(
      (Success("success key"), Failure(new RuntimeException("failed value"))), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 1
    exception.throwables.head.getMessage should be("failed value")
  }

  it should "sequence a failed key with a successful value" in {
    val result: Try[(String, String)] = sequenceTuple(
      (Failure(new RuntimeException("failed key")), Success("success value")), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 1
    exception.throwables.head.getMessage should be("failed key")
  }

  it should "sequence a failed key with a failed value" in {
    val result: Try[(String, String)] = sequenceTuple(
      (Failure(new RuntimeException("failed key")), Failure(new RuntimeException("failed value"))), "prefix")
    result.isFailure should be(true)
    result.failed.get should be(an[AggregatedException])
    val exception = result.failed.get.asInstanceOf[AggregatedException]
    exception.exceptionContext should be("prefix")
    exception.throwables should have size 2
    exception.throwables.head.getMessage should be("failed key")
    exception.throwables.tail.head.getMessage should be("failed value")
  }
}
