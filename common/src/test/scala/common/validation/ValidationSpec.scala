package common.validation

import cats.data.Validated.{Invalid, Valid}
import cats.data.{NonEmptyList, ValidatedNel}
import cats.syntax.validated._
import common.exception.AggregatedMessageException
import common.validation.Validation._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.slf4j.Logger
import org.specs2.mock.Mockito

import scala.util.{Failure, Success}


class ValidationSpec extends AnyFlatSpec with Matchers with Mockito {

  behavior of "Validation"

  it should "warn unrecognized keys" in {
    var warnings = List.empty[Any]
    val mockLogger = mock[Logger]
    mockLogger.warn(anyString).answers(warnings :+= _)
    val keys = Set("hello")
    val reference = Set("world")
    val context = "warnings"
    warnNotRecognized(keys, reference, context, mockLogger)
    warnings should contain theSameElementsAs List("Unrecognized configuration key(s) for warnings: hello")
  }

  it should "not warn recognized keys" in {
    var warnings = List.empty[Any]
    val mockLogger = mock[Logger]
    mockLogger.warn(anyString).answers(warnings :+= _)
    val keys = Set("hello")
    val reference = Set("hello", "world")
    val context = "warnings"
    warnNotRecognized(keys, reference, context, mockLogger)
    warnings should be(empty)
  }

  it should "convert a valid to a successful future" in {
    val nel: ValidatedNel[String, String] = "ok".valid
    val actual = nel.toFuture(_ => throw new RuntimeException("Should not be used"))
    actual.value should be(Option(Success("ok")))
  }

  it should "convert a invalidNel to an unsuccessful future" in {
    val nel: ValidatedNel[String, String] = "failed".invalidNel
    val actual = nel.toFuture(nel => new RuntimeException(nel.toList.mkString("processed ", ",", " and caught")))
    actual.value match {
      case Some(Failure(exception: RuntimeException)) => exception.getMessage should be("processed failed and caught")
      case other => fail(s"expected message: 'processed failed and caught' but got '$other'")
    }
  }

  it should "succeed to validate a valid value" in {
    val result = validate("hello")
    result should be("hello".valid)
  }

  it should "fail to validate an invalid value" in {
    val result = validate(throw new RuntimeException("fail"))
    result should be("fail".invalidNel)
  }

  it should "convert a Try to an ErrorOr" in {
    val success = Success("yeah")
    val failure = Failure(new Exception(":("))
    import common.validation.Validation._
    success.toErrorOr shouldBe Valid("yeah")
    failure.toErrorOr shouldBe Invalid(NonEmptyList.of(":("))
  }

  it should "convert a Try to a Checked" in {
    val success = Success("yeah")
    val failure = Failure(new Exception(":("))
    import common.validation.Validation._
    success.toChecked shouldBe Right("yeah")
    failure.toChecked shouldBe Left(NonEmptyList.of(":("))
  }

  it should "convert an ErrorOr to a Try" in {
    val valid = "yeah".valid
    val invalid = ":(".invalidNel
    import common.validation.Validation._
    valid.toTry should be(Success("yeah"))
    val exception = intercept[AggregatedMessageException](invalid.toTry.get)
    exception.exceptionContext should be("Error(s)")
    exception.errorMessages should contain theSameElementsAs List(":(")
  }

  it should "convert a Some to an ErrorOr" in {
    Option("ok").toErrorOr("not used") should be("ok".valid)
  }

  it should "convert a None to an ErrorOr" in {
    None.toErrorOr("error message") should be("error message".invalidNel)
  }

  it should "convert a Some to an Checked" in {
    Option("ok").toChecked("not used") should be(Right("ok"))
  }

  it should "convert a None to an Checked" in {
    None.toChecked("error message") should be(Left(NonEmptyList.one("error message")))
  }

}
