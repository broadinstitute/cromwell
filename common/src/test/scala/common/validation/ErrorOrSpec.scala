package common.validation

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import common.assertion.CromwellTimeoutSpec
import common.validation.ErrorOr._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class ErrorOrSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

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

  it should "flatMapN 1-tuples into a Valid string" in {
    val result = Tuple1(valid0).flatMapN(Array(_).mkString.valid)
    result should be(Valid("0"))
  }

  it should "flatMapN 2-tuples into a Valid string" in {
    val result = (valid0, valid1).flatMapN(Array(_, _).mkString.valid)
    result should be(Valid("01"))
  }

  it should "flatMapN 3-tuples into a Valid string" in {
    val result = (valid0, valid1, valid2).flatMapN(Array(_, _, _).mkString.valid)
    result should be(Valid("012"))
  }

  it should "flatMapN 4-tuples into a Valid string" in {
    val result = (valid0, valid1, valid2, valid0).flatMapN(Array(_, _, _, _).mkString.valid)
    result should be(Valid("0120"))
  }

  it should "flatMapN 5-tuples into a Valid string" in {
    val result = (valid0, valid1, valid2, valid0, valid1).flatMapN(Array(_, _, _, _, _).mkString.valid)
    result should be(Valid("01201"))
  }

  it should "flatMapN 6-tuples into a Valid string" in {
    val result = (valid0, valid1, valid2, valid0, valid1, valid2).flatMapN(Array(_, _, _, _, _, _).mkString.valid)
    result should be(Valid("012012"))
  }

  it should "flatMapN 7-tuples into a Valid string" in {
    val result = (valid0, valid1, valid2, valid0, valid1, valid2, valid0)
      .flatMapN(Array(_, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("0120120"))
  }

  it should "flatMapN 8-tuples into a Valid string" in {
    val result = (valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1)
      .flatMapN(Array(_, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("01201201"))
  }

  it should "flatMapN 9-tuples into a Valid string" in {
    val result = (valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("012012012"))
  }

  it should "flatMapN 10-tuples into a Valid string" in {
    val result = (valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("0120120120"))
  }

  it should "flatMapN 11-tuples into a Valid string" in {
    val result = (valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("01201201201"))
  }

  it should "flatMapN 12-tuples into a Valid string" in {
    val result = (valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("012012012012"))
  }

  it should "flatMapN 13-tuples into a Valid string" in {
    val result = (
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2,
      valid0)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("0120120120120"))
  }

  it should "flatMapN 14-tuples into a Valid string" in {
    val result = (
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2,
      valid0, valid1)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("01201201201201"))
  }

  it should "flatMapN 15-tuples into a Valid string" in {
    val result = (
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2,
      valid0, valid1, valid2)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("012012012012012"))
  }

  it should "flatMapN 16-tuples into a Valid string" in {
    val result = (
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2,
      valid0, valid1, valid2, valid0)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("0120120120120120"))
  }

  it should "flatMapN 17-tuples into a Valid string" in {
    val result = (
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2,
      valid0, valid1, valid2, valid0, valid1)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("01201201201201201"))
  }

  it should "flatMapN 18-tuples into a Valid string" in {
    val result = (
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2,
      valid0, valid1, valid2, valid0, valid1, valid2)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("012012012012012012"))
  }

  it should "flatMapN 19-tuples into a Valid string" in {
    val result = (
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2,
      valid0, valid1, valid2, valid0, valid1, valid2, valid0)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("0120120120120120120"))
  }

  it should "flatMapN 20-tuples into a Valid string" in {
    val result = (
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2,
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("01201201201201201201"))
  }

  it should "flatMapN 21-tuples into a Valid string" in {
    val result = (
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2,
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2)
      .flatMapN(Array(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("012012012012012012012"))
  }

  it should "flatMapN 22-tuples into a Valid string" in {
    val result = (
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2,
      valid0, valid1, valid2, valid0, valid1, valid2, valid0, valid1, valid2, valid0).flatMapN(
      Array(_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _).mkString.valid)
    result should be(Valid("0120120120120120120120"))
  }

  it should "short-circuit when given Invalid inputs to flatMapN" in {
    (((valid2, valid0) flatMapN errorOrDiv, valid1) flatMapN errorOrDiv) should be(DivBy0Error.invalidNel)
  }

  it should "generate a flatMapN" in {
    val result = ErrorOrGen.mkShortCircuitingFlatMapTupleNFunction(1) should be(
      """|implicit class ShortCircuitingFlatMapTuple1[A](val t1: (ErrorOr[A])) extends AnyVal {
         |    def flatMapN[T_OUT](f1: (A) => ErrorOr[T_OUT]): ErrorOr[T_OUT] = t1.tupled flatMap f1.tupled
         |}
         |
         |""".stripMargin)
    result
  }

}
