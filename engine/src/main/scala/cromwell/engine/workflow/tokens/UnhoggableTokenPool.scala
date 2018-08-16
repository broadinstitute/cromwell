package cromwell.engine.workflow.tokens

import java.util.UUID
import java.util.concurrent.atomic.AtomicBoolean

import cromwell.core.JobExecutionToken
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import cromwell.engine.workflow.tokens.UnhoggableTokenPool.{ComeBackLater, Oink, TokenHoggingLease, UnhoggableTokenPoolResult}
import io.github.andrebeat.pool._

import scala.collection.immutable.HashSet
import scala.collection.mutable

final case class UnhoggableTokenPool(tokenType: JobExecutionTokenType) extends SimplePool[JobExecutionToken](
  capacity = tokenType.maxPoolSize.getOrElse(UnhoggableTokenPool.MaxCapacity),
  referenceType = ReferenceType.Strong,
  _factory = () => JobExecutionToken(tokenType, UUID.randomUUID()),
  _reset = Function.const(()),
  _dispose = Function.const(()),
  _healthCheck = Function.const(true)) {

  lazy val hogLimitOption = tokenType.maxPoolSize map { limit => math.round(limit.floatValue() / tokenType.hogFactor.floatValue()) }

  private[this] val hogGroupAssignments: mutable.Map[String, HashSet[JobExecutionToken]] = mutable.Map.empty
  private[this] var hogGroupAssignmentCount = 0

  override def tryAcquire(): Option[Lease[JobExecutionToken]] = throw new UnsupportedOperationException("Use tryAcquire(hogGroup)")

  def tryAcquire(hogGroup: String): UnhoggableTokenPoolResult = {

    hogLimitOption match {
      case Some(hogLimit) =>
        if (hogGroupAssignmentCount + 1 > hogLimit) {
          Oink
        } else {
          synchronized {
            val thisHogSet = hogGroupAssignments.getOrElse(hogGroup, HashSet.empty)
            if (thisHogSet.size + 1 < hogLimit) {
              super.tryAcquire() match {
                case Some(lease) =>
                  val hoggingLease = TokenHoggingLease(lease, this)
                  hogGroupAssignments += (hogGroup -> thisHogSet.+(hoggingLease.get))
                  hoggingLease
                case None => ComeBackLater
              }
            } else {
              Oink
            }
          }
        }
      case None =>
        super.tryAcquire() match {
          case Some(lease) =>
            val hoggingLease = TokenHoggingLease(lease, this)
            hoggingLease
          case None => ComeBackLater
        }
    }
  }

  def unhog(hogGroup: String, lease: Lease[JobExecutionToken]): Unit = {
    hogLimitOption foreach { _ =>
      synchronized {
        hogGroupAssignments += (hogGroup -> hogGroupAssignments.getOrElse(hogGroup, HashSet.empty).-(lease.get))
      }
    }
  }
}

object UnhoggableTokenPool {
  // 10 Million
  val MaxCapacity = 10 * 1000 * 1000

  trait UnhoggableTokenPoolResult

  final case class TokenHoggingLease(lease: Lease[JobExecutionToken], pool: UnhoggableTokenPool) extends Lease[JobExecutionToken] with UnhoggableTokenPoolResult {
    private[this] val dirty = new AtomicBoolean(false)
    override protected[this] def a: JobExecutionToken = lease.get

    override protected[this] def handleRelease(): Unit = {
      if (dirty.compareAndSet(false, true)) {

      }
      lease.release()
    }
    override protected[this] def handleInvalidate(): Unit = lease.invalidate()
  }

  /**
    * You didn't get a lease because the pool is empty
    */
  case object ComeBackLater extends UnhoggableTokenPoolResult

  /**
    * You didn't get a lease... because you're being a hog
    */
  case object Oink extends UnhoggableTokenPoolResult

}
