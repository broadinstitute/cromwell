package cromwell.core.io

import scala.util.{Failure, Success, Try}

/**
  * Generic trait for values returned after a command is executed. Can be Success or Failure.
  *
  * @tparam T type of the returned value if success
  */
sealed trait IoAck[T] {
  /**
    * Original command
    */
  def command: IoCommand[T]
  def toTry: Try[T]
}

case class IoSuccess[T](command: IoCommand[T], result: T) extends IoAck[T] {
  override def toTry = Success(result)
}

trait IoFailAck[T] extends IoAck[T] {
  val failure: Throwable
  override def toTry = Failure(failure)
}

/** Failure of an unspecified variety. */
case class IoFailure[T](command: IoCommand[T], override val failure: Throwable) extends IoFailAck[T]
/** Failure of the 403 persuasion. */
case class IoForbiddenFailure[T](command: IoCommand[T], override val failure: Throwable, forbiddenPath: String) extends IoFailAck[T]
