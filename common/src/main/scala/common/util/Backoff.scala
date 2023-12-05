package common.util

import scala.concurrent.duration.FiniteDuration

trait Backoff {

  /** Next interval in millis */
  def backoffMillis: Long

  /** Get the next instance of backoff. This should be called after every call to backoffMillis */
  def next: Backoff
}

object Backoff {
  def staticBackoff(time: FiniteDuration) = StaticBackoff(time)
}

case class StaticBackoff(time: FiniteDuration) extends Backoff {
  override def backoffMillis = time.toMillis
  override def next = this
}
