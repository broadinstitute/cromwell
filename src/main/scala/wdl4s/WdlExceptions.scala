package wdl4s

import scalaz.NonEmptyList

trait ThrowableWithErrors extends Exception {
  def message: String
  def errors: NonEmptyList[String]

  override def getMessage = {
    s"""$message\n${errors.list.mkString("\n")}"""
  }
}

class UnsatisfiedInputsException(message: String) extends Exception(message)
class ValidationException(val message: String, val errors: NonEmptyList[String]) extends ThrowableWithErrors
