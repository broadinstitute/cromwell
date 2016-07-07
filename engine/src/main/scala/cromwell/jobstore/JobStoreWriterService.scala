package cromwell.jobstore

import akka.actor.{Actor, ActorRef}
import com.typesafe.config.Config
import cromwell.core.WorkflowId
import cromwell.jobstore.JobStoreWriterService.JobStoreWriterServiceCommand
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

import scala.concurrent.Future

/**
  * Joins the service registry API to the JobStoreWriterActor.
  *
  * This level of indirection is a tiny bit awkward but allows the database to be injected.
  */
case class JobStoreWriterService(serviceConfig: Config, globalConfig: Config) extends Actor {

  // TODO: Replace with a real database, probably from config.
  val database = new ConsoleOutputJobStoreDatabase()
  var JSWActor = context.actorOf(JobStoreWriterActor.props(database))

  override def receive: Receive = {
    case command: JobStoreWriterServiceCommand => JSWActor.tell(command, sender)
  }
}

object JobStoreWriterService {
  sealed trait JobStoreWriterServiceCommand extends ServiceRegistryMessage { override def serviceName = "JobStoreWriter" }
  case class RegisterJobCompleted(jobKey: JobStoreKey, jobResult: JobResult) extends JobStoreWriterServiceCommand
  case class RegisterWorkflowCompleted(workflowId: WorkflowId) extends JobStoreWriterServiceCommand

  sealed trait JobStoreWriterServiceResponse
  case class JobStoreWriteSuccess(originalCommand: JobStoreWriterServiceCommand) extends JobStoreWriterServiceResponse
  case class JobStoreWriteFailure(originalCommand: JobStoreWriterServiceCommand, reason: Throwable) extends JobStoreWriterServiceResponse
}

class ConsoleOutputJobStoreDatabase extends JobStoreDatabase {

  override def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId]): Future[Unit] = {
    System.out.println(s"Jobs completed: ${jobCompletions.mkString}")
    System.out.println(s"Workflows completed: ${workflowCompletions.mkString}")
    Future.successful(())
  }
}