package cromwell.backend.impl.jes.statuspolling

import akka.actor.{Actor, ActorLogging, ActorRef}
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.RunPipelineRequest
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.{JesApiRunCreationQueryFailed}
import cromwell.backend.standard.StandardAsyncJob

import scala.concurrent.{Future, Promise}
import scala.util.{Failure, Success, Try}

/**
  * I'm putting this stuff in a mixin to avoid polluting the main class.
  *
  * Be sure to make the main class's receive look like:
  *   override def receive = runCreationClientReceive orElse { ... }
  */
trait JesRunCreationClient { this: Actor with ActorLogging =>

  private var runCreationClientPromise: Option[Promise[StandardAsyncJob]] = None

  val pollingActor: ActorRef

  def runCreationClientReceive: Actor.Receive = {
    case job: StandardAsyncJob =>
      completePromise(Success(job))
    case JesApiRunCreationQueryFailed(_, e) => completePromise(Failure(e))
  }

  private def completePromise(job: Try[StandardAsyncJob]) = {
    runCreationClientPromise foreach { _.complete(job) }
    runCreationClientPromise = None
  }

  def runPipeline(genomicsInterface: Genomics, rpr: RunPipelineRequest): Future[StandardAsyncJob] = {
    runCreationClientPromise match {
      case Some(p) =>
        p.future
      case None =>
        pollingActor ! JesApiQueryManager.DoCreateRun(genomicsInterface, rpr)
        val newPromise = Promise[StandardAsyncJob]()
        runCreationClientPromise = Option(newPromise)
        newPromise.future
    }
  }
}
