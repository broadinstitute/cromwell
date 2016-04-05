package cromwell.services

import akka.actor.Actor
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowMetadata
import cromwell.engine.{ExecutionEventEntry, WorkflowId}
import cromwell.engine.backend.JobKey
import cromwell.services.MetadataServiceActor._
import cromwell.webservice.WorkflowMetadataResponse
import org.joda.time.DateTime

object MetadataServiceActor {

  sealed trait MetadataServiceActorMessage
  case class SetStartTime(jobKey: JobKey, startTime: DateTime)
  case class SetEndTime(jobKey: JobKey, time: DateTime)
  case class AddExecutionEvent(jobKey: JobKey, executionEventEntry: ExecutionEventEntry)
  /*
  We'd have a "set" message for each of these:

  Workflow Level:
                        id: String,
                        workflowName: String,
                        status: String,
                        submission: DateTime,
                        start: Option[DateTime],
                        end: Option[DateTime],
                        inputs: JsObject,
                        outputs: Option[Map[String, WdlValue]],
                        failures: Option[Seq[FailureEventEntry]]
  Call Level:
                        inputs: Map[String, WdlValue],
                        executionStatus: String,
                        backend: Option[String],
                        backendStatus: Option[String],
                        outputs: Option[Map[String, WdlValue]],
                        start: Option[DateTime],
                        end: Option[DateTime],
                        jobId: Option[String],
                        returnCode: Option[Int],
                        shardIndex: Int,
                        stdout: Option[WdlFile],
                        stderr: Option[WdlFile],
                        backendLogs: Option[Map[String, WdlFile]],
                        executionEvents: Seq[ExecutionEventEntry],
                        attempt: Int,
                        runtimeAttributes: Map[String, String],
                        preemptible: Option[Boolean],
                        failures: Option[Seq[FailureEventEntry]])
   */

}

class MetadataServiceActor extends Actor {

  def receive = {
    case SetStartTime(jobKey, startTime) => ??? // Write this start time information to the database

      // This replaces the functionality on WorkflowManagerActor - the slick endpoint would use this service actor instead.
      // NB: Move the messages WorkflowMetadata and WorkflowMetadataResponse into the MetadataServiceActor object instead.
    case WorkflowMetadata(workflowId) =>
      val workflowMetadataResponse: WorkflowMetadataResponse = ??? // Build it from the database
      sender ! workflowMetadataResponse
  }

}
