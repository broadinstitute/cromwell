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
case class IoFailure[T](command: IoCommand[T], failure: Throwable) extends IoAck[T] {
  override def toTry = Failure(failure)
}

case class IoRetry[T](command: IoCommand[T], failure: Throwable) extends IoAck[T] {
  override def toTry = Failure(failure)
}
