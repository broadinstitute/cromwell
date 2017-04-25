package cromwell.core.actor

import cats.data.NonEmptyVector
import org.slf4j.LoggerFactory

import scala.util.{Failure, Success, Try}


/** A collection of state, data, and message types to support batched database writes. */
object BatchingDbWriter {
  val logger = LoggerFactory.getLogger("BatchingDbWriteActor")

  /** Data for batched database writes. */
  sealed trait BatchingDbWriterData {
    def addData[D](datum: D): BatchingDbWriterData = addData(Vector(datum))
    def addData[D](data: Iterable[D]): BatchingDbWriterData = {
      Try(NonEmptyVector.fromVector(data.toVector)) match {
        case Success(Some(v)) =>
          val newEvents = this match {
            case NoData => v
            case HasData(e) => e.concatNev(v)
          }
          HasData(newEvents)
        case Success(None) => this
        case Failure(f) =>
          val dataSample = data.take(3).mkString(", ") + (if (data.size > 3) ", ..." else "")
          logger.error(s"Failed processing batched data: $dataSample. Data will be dropped and not be sent to the database.", f)
          this
      }
    }

    def length: Int = this match {
      case NoData => 0
      case HasData(e) => e.length
    }
  }

  case object NoData extends BatchingDbWriterData
  case class HasData[E](events: NonEmptyVector[E]) extends BatchingDbWriterData

  /** The states for batched database writes. */
  sealed trait BatchingDbWriterState
  case object WaitingToWrite extends BatchingDbWriterState
  case object WritingToDb extends BatchingDbWriterState

  /** The message types for batched database writes. */
  sealed trait BatchingDbWriterMessage
  case object DbWriteComplete extends BatchingDbWriterMessage
  case object FlushBatchToDb extends BatchingDbWriterMessage
  case object ScheduledFlushToDb extends BatchingDbWriterMessage
}
