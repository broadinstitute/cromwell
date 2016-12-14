package lenthall.validation

import cats.syntax.validated._
import lenthall.validation.ErrorOr._
import org.scalatest.{FlatSpec, Matchers}

class ErrorOrSpec extends FlatSpec with Matchers {

  behavior of "ErrorOr"

  it should "flatMap valid nels" in {
    val errorOrA: ErrorOr[String] = "hello".valid
    val errorOrB: ErrorOr[String] = "world".valid

    val errorMapped = new ShortCircuitingFlatMap(errorOrA).flatMap(_ => errorOrB)
    errorMapped should be("world".valid)
  }

  it should "flatMap invalid nels" in {
    val errorOrA: ErrorOr[String] = "hello".invalidNel
    val errorOrB: ErrorOr[String] = "world".valid

    val errorMapped = new ShortCircuitingFlatMap(errorOrA).flatMap(_ => errorOrB)
    errorMapped should be("hello".invalidNel)
  }

}
