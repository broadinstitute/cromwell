package cromwell.core.actor

import cats.data.NonEmptyVector

import cats.syntax.traverse._
import cats.instances.vector._
import cats.instances.future._
import scala.concurrent.Future
import scala.concurrent.duration.Duration

/**
  * By removing the periodic flush and setting the batch size to 1,
  * this actor modifies the behavior of BatchActor so that it throttles processing commands one at a time.
  */
abstract class ThrottlerActor[C] extends BatchActor[C](Duration.Zero, 1) {
  override def weightFunction(command: C) = 1
  override final def process(data: NonEmptyVector[C]): Future[Int] = {
    // This ShouldNotBePossible™ but in case it happens, instead of dropping elements process them all anyway.
    // Explanation: batch size is 1 which means as soon as we receive 1 element, the process method should be called.
    // Because the BatchActor calls the process method with vector of elements which total weight is batch size, and because
    // each element's weight is 1, this method should always be called with a vector of size 1. If that is not the case it means
    // there's a bug in the batch actor.
    if (data.tail.nonEmpty) {
      log.error("{} is throttled and is not supposed to process more than one element at a time !", self.path.name)
      data.toVector.traverse[Future, Any](processHead).map(_.length)
    } else processHead(data.head).map(_ => 1)
  }
  def processHead(head: C): Future[Any]
}
