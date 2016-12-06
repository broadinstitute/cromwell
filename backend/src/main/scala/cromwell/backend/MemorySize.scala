package cromwell.backend


import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._
import lenthall.validation.ErrorOr._
import mouse.string._

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import wdl4s.parser.MemoryUnit


object MemorySize {
  val memoryPattern = """(\d+(?:\.\d+)?)\s*(\w+)""".r

  def parse(unparsed: String): Try[MemorySize] = {
    unparsed match {
      case memoryPattern(amountString, unitString) =>
        val amount: ErrorOr[Double] = amountString.parseDouble leftMap {
          _.getMessage
        } toValidatedNel
        val unit: ErrorOr[MemoryUnit] = MemoryUnit.values find {
          _.suffixes.contains(unitString)
        } match {
          case Some(s) => s.validNel
          case None => s"$unitString is an invalid memory unit".invalidNel
        }
        (amount |@| unit) map { (a, u) => new MemorySize(a, u) } match {
          case Valid(memorySize) => Success(memorySize)
          case Invalid(nel) => Failure(new UnsupportedOperationException(nel.toList.mkString("\n")))
        }
      case _ => Failure(new UnsupportedOperationException(s"$unparsed should be of the form 'X Unit' where X is a number, e.g. 8 GB"))
    }
  }
}

case class MemorySize(amount: Double, unit: MemoryUnit) {
  def bytes: Double = amount * unit.bytes

  def to(unit: MemoryUnit): MemorySize = MemorySize(this.bytes / unit.bytes, unit)

  override def toString: String = {
    val adjustedAmount = (unit, amount) match {
      case (MemoryUnit.Bytes, a) => a.ceil.toLong.toString
      case (_, a) if a == a.toLong => a.toLong.toString
      case (_, a) => a.toString
    }
    s"$adjustedAmount ${unit.suffixes(0)}"
  }
}
