package cromwell.backend.impl.aws

import java.util.concurrent.Executors

import akka.actor.{Actor, ActorLogging, Props}
import cats.effect.{IO, Timer}
import cromwell.backend.impl.aws.GoSlowJobSubmitActor.{CheckForWork, SubmitForMe}
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.batch.model.SubmitJobResponse

import scala.concurrent.{ExecutionContext, Promise}
import scala.concurrent.duration._
import scala.util.{Failure, Success}

class GoSlowJobSubmitActor(configRegion: Option[Region]) extends Actor with ActorLogging {

  implicit val ec: ExecutionContext = context.dispatcher
  implicit val timer: Timer[IO] = cats.effect.IO.timer(ExecutionContext.fromExecutor(Executors.newSingleThreadScheduledExecutor))

  val WorkInterval = 100.millis

  lazy val client = {
    val builder = BatchClient.builder()
    configRegion.foreach(builder.region)
    builder.build
  }

  // Maps job ID to status
  var workQueue: Vector[SubmitForMe] = Vector.empty


  override def receive = {
    case sfm: SubmitForMe =>
      workQueue :+= sfm

    case CheckForWork => workQueue.headOption match {
      case Some(SubmitForMe(batchJob, attributes, completionPromise)) =>
        batchJob.submitJob[IO]().run(attributes).unsafeToFuture() onComplete {
          case Success(value) =>
            completionPromise.success(value)
            scheduleWorkCheck(WorkInterval)
          case Failure(error) =>
            completionPromise.failure(error)
            log.error(error, s"Submission of new AWS job failed")
            scheduleWorkCheck(WorkInterval)
        }
        workQueue = workQueue.tail
      case None =>
        scheduleWorkCheck(WorkInterval)
    }
  }

  def scheduleWorkCheck(in: FiniteDuration): Unit = {
    context.system.scheduler.scheduleOnce(in) { self ! CheckForWork }
    ()
  }

  scheduleWorkCheck(WorkInterval)
}

object GoSlowJobSubmitActor {

  sealed trait GoSlowSubmitActorMessage

  case object CheckForWork

  final case class SubmitForMe(job: AwsBatchJob, attributes: AwsBatchAttributes, completionPromise: Promise[SubmitJobResponse]) extends GoSlowSubmitActorMessage

  def props(configRegion: Option[Region]) = Props(new GoSlowJobSubmitActor(configRegion))
}
