package wdl4s

import cats.data.NonEmptyList

trait ExceptionWithErrors extends Exception {
  def message: String
  def errors: NonEmptyList[String]

  override def getMessage = {
    s"""$message\n${errors.toList.mkString("\n")}"""
  }
}

class UnsatisfiedInputsException(message: String) extends Exception(message)
class ValidationException(val message: String, val errors: NonEmptyList[String]) extends ExceptionWithErrors
