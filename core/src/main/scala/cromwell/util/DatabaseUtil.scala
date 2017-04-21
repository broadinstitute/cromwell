package cromwell.util

import java.sql.SQLTransientException

import akka.NotUsed
import akka.actor.ActorSystem
import akka.stream.scaladsl.Source
import cromwell.core.retry.{Retry, SimpleExponentialBackoff}

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._
import scala.language.postfixOps

object DatabaseUtil {
  private def isTransient(throwable: Throwable): Boolean = throwable match {
    case _: SQLTransientException => true
    case _ => false
  }

  def withRetry[A](f: () => Future[A])(implicit actorSystem: ActorSystem): Future[A] = {
    val RetryBackoff = SimpleExponentialBackoff(50 millis, 1 seconds, 1D)
    Retry.withRetry(f, maxRetries = Option(10), backoff = RetryBackoff, isTransient = isTransient)
  }

  def paginatedSource[T](pageSize: Int, queryFunction: (Int, Int) => Future[Seq[T]])(implicit ec: ExecutionContext): Source[Seq[T], NotUsed] = {
    Source.unfoldAsync((1, 0))({
      case (pageNumber, previousBatchSize) if pageNumber == 1 || previousBatchSize == pageSize => 
        queryFunction(pageSize, pageNumber).map(contents => Some((pageNumber + 1, contents.length) -> contents))
      case _ => Future.successful(None)
    })
  }

  def oneByOneSource[T](queryFunction: (Int) => Future[Option[T]])(implicit ec: ExecutionContext): Source[T, NotUsed] = {
    Source.unfoldAsync(1)({ pageNumber =>
      queryFunction(pageNumber) map {
        case Some(next) => Option(pageNumber + 1 -> next)
        case None => None
      }
    })
  }

}
