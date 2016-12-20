package cromwell.services.io.message

import cromwell.services.io.IoActorCommand

import scala.util.{Failure, Success, Try}

/**
  * Generic trait for values returned after a command is executed. Can be Success or Failure.
  *
  * @tparam T type of the returned value if success
  */
sealed trait IoActorMessageAck[T] {
  /**
    * Original command
    */
  def command: IoActorCommand[T]

  /**
    * @return Success(value) in case of success or Failure(exception) in case of failure
    */
  def toOption: Try[T]
}

case class IoSuccess[T](command: IoActorCommand[T], result: T) extends IoActorMessageAck[T] {
  override def toOption = Success(result)
}

case class IoRetried[T](command: IoActorCommand[T], failure: Throwable) extends IoActorMessageAck[T] {
  override def toOption = Failure(failure)
}

case class IoFailure[T](command: IoActorCommand[T], failure: Throwable) extends IoActorMessageAck[T] {
  override def toOption = Failure(failure)
}
