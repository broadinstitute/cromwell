package cromwell.database

import akka.actor.ActorSystem
import cromwell.core.retry.{Retry, SimpleExponentialBackoff}
import cromwell.database.slick.SlickDatabase

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps


object CromwellDatabase {
  val databaseInterface: SqlDatabase = new SlickDatabase()
}

trait Database {
  def databaseInterface: SqlDatabase

  def withRetry[A](f: => Future[A])(implicit actorSystem: ActorSystem): Future[A] = {
    val RetryBackoff = SimpleExponentialBackoff(50 millis, 1 seconds, 1D)
    Retry.withRetry(() => f, maxRetries = Option(10), backoff = RetryBackoff, isTransient = databaseInterface.isTransient)
  }
}

trait CromwellDatabase extends Database {
  override def databaseInterface: SqlDatabase = CromwellDatabase.databaseInterface
}
