package lenthall.validation

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.data.Validated.Valid
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

  it should "sequence a map of valid ErrorOrs" in {
    val eo1 = Valid("good1")
    val eo2 = Valid("good2")
    Map(1 -> eo1, 2 -> eo2).sequence should be(Valid(Map(1 -> "good1", 2 -> "good2")))
  }

  it should "sequence a map of mixed ErrorOrs" in {
    val eo1 = Valid("good1")
    val eo2 = "bad2".invalidNel
    val eo3 = Valid("good3")
    val eo4 = "bad4".invalidNel
    Map(1 -> eo1, 2 -> eo2, 3 -> eo3, 4 -> eo4).sequence should be(Invalid(NonEmptyList("bad2", List("bad4"))))
  }

  val DivBy0Error: String = "Divide by 0!"
  def errorOrDiv(v1: Int, v2: Int): ErrorOr[Double] = if (v2 != 0) { Valid(v1.toDouble / v2.toDouble) } else { DivBy0Error.invalidNel }
  def errorOrDiv(v1: Double, v2: Int): ErrorOr[Double] = if (v2 != 0) { Valid(v1.toDouble / v2.toDouble) } else { DivBy0Error.invalidNel }
  def errorOrSelect(v1: Int, v2: Int, v3: Int, v4: Int, v5: Int, v6: Int, v7: Int,
                    v8: Int, v9: Int, v10: Int, v11: Int, v12: Int, v13: Int, v14: Int,
                    v15: Int, v16: Int, v17: Int, v18: Int, v19: Int, v20: Int, v21: Int, v22: Int): ErrorOr[Int] = Valid(v4 + v6 + v22)

  val valid0 = Valid(0)
  val valid1 = Valid(1)
  val valid2 = Valid(2)

  it should "flatMapN 2-tuples of Valid's into a Valid" in {
    (valid2, valid1) flatMapN { (x, y) => errorOrDiv(x, y) } should be(Valid(2d))
  }

  it should "flatMapN 2-tuples of Valid's into an invalid" in {
    (valid2, valid0) flatMapN errorOrDiv should be(DivBy0Error.invalidNel)
  }

  it should "flatMapN 22-tuples into a Valid" in {
    (valid0, valid1, valid2,
      valid0, valid1, valid2,
      valid0, valid1, valid2,
      valid0, valid1, valid2,
      valid0, valid1, valid2,
      valid0, valid1, valid2,
      valid0, valid1, valid2, valid0) flatMapN errorOrSelect should be(Valid(0 + 2 + 0))
  }

  it should "short-circuit when given Invalid inputs to flatMapN" in {
    (((valid2, valid0) flatMapN errorOrDiv, valid1) flatMapN errorOrDiv) should be(DivBy0Error.invalidNel)
  }

}
