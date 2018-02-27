package cromwell.core.actor

import cats.data.NonEmptyVector

import cats.syntax.traverse._
import cats.instances.vector._
import cats.instances.future._
import scala.concurrent.Future
import scala.concurrent.duration.Duration

/**
  * By removing the periodic flush and setting the batch size to 1,
  * This transforms the BatchActor into throttler processing commands one at a time.
  */
abstract class ThrottlerActor[C] extends BatchActor[C](Duration.Zero, 1) {
  override def weightFunction(command: C) = 1
  override final def process(data: NonEmptyVector[C]): Future[Int] = {
    // This ShouldNotBePossible™ but in case it happens, instead of dropping elements process them all anyway
    if (data.tail.nonEmpty) {
      log.error("{} is a read actor and is not supposed to process more than one element at a time !", self.path.name)
      data.toVector.traverse[Future, Any](processHead).map(_.length)
    } else processHead(data.head).map(_ => 1)
  }
  def processHead(head: C): Future[Any]
}
