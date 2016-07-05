package cromwell.jobstore

import akka.actor.{Actor, ActorLogging, Props}
import akka.event.LoggingReceive
import cromwell.backend.BackendJobDescriptorKey

object JobStoreReader {
  def props() = Props(new JobStoreReader())
}

class JobStoreReader extends Actor with ActorLogging {
  override def receive = LoggingReceive {
    case QueryJobCompletion(jobKey) => sender ! JobNotComplete(jobKey)
    case unknownMessage => log.error(s"Unexpected message to JobStoreReader: $unknownMessage")
  }
}

/**
  * Message to query the JobStoreReader, asks whether the specified job has already been completed.
  */
case class QueryJobCompletion(jobKey: BackendJobDescriptorKey)

/**
  * Message which indicates that a job has already completed, and contains the results of the job
  */
case class JobComplete(jobKey: BackendJobDescriptorKey, jobResult: JobResult)

/**
  * Indicates that the job has not been completed yet. Makes no statement about whether the job is
  * running versus unstarted or (maybe?) doesn't even exist!
  */
case class JobNotComplete(jobKey: BackendJobDescriptorKey)
