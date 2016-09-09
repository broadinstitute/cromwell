package cromwell.backend

import wdl4s.parser.MemoryUnit

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz.Scalaz._

object MemorySize {
  val memoryPattern = """(\d+(?:\.\d+)?)\s*(\w+)""".r

  def parse(unparsed: String): Try[MemorySize] = {
    unparsed match {
      case memoryPattern(amountString, unitString) =>
        val amount = amountString.parseDouble leftMap {
          _.getMessage
        } toValidationNel
        val unit = MemoryUnit.values find {
          _.suffixes.contains(unitString)
        } match {
          case Some(s) => s.successNel[String]
          case None => s"$unitString is an invalid memory unit".failureNel
        }
        (amount |@| unit) { (a, u) => new MemorySize(a, u) } match {
          case scalaz.Success(memorySize) => Success(memorySize)
          case scalaz.Failure(nel) => Failure(new UnsupportedOperationException(nel.list.toList.mkString("\n")))
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
