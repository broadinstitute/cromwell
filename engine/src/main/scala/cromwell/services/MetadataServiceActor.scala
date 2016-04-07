package cromwell.services

import akka.actor.Actor
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowMetadata
import cromwell.engine.{ExecutionEventEntry, WorkflowId}
import cromwell.engine.backend.JobKey
import cromwell.services.MetadataServiceActor._
import cromwell.webservice.WorkflowMetadataResponse
import org.joda.time.DateTime
import wdl4s.values.WdlValue

object MetadataServiceActor {

  sealed trait MetadataEvent[V] {
    def workflowId: WorkflowId
    def jobKey: JobKey
    def eventName: String
    def eventTime: DateTime
    def value: V
  }

  /** Every instance of this received by the MetadataServiceActor is simply appended to the DB table; previously existing rows
    * for this workflowId + jobKey are not modified.  This retains the full record of what changed and when.  `submission`, `start`, and `end` might
    * be derived from the appropriate rows with eventName = "Workflow Status".  */
  final case class WorkflowStatusEvent(workflowId: WorkflowId, jobKey: JobKey, eventTime: DateTime, value: String) extends MetadataEvent[String] {
    override def eventName = "Workflow Status"
  }

  /**
    * Note the additional fqn argument in this one:
    */
  final case class WorkflowOutputGenerated(workflowId: WorkflowId, jobKey: JobKey, eventTime: DateTime, fqn: String, value: WdlValue) extends MetadataEvent[WdlValue] {
    override def eventName = "Output"
  }

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
    case WorkflowStatusEvent(workflowId: WorkflowId, jobKey: JobKey, eventTime: DateTime, value: String) =>
      // Append the status change event to the DB's event journal (wrong word?)
      ???
    case WorkflowOutputGenerated(workflowId: WorkflowId, jobKey: JobKey, eventTime: DateTime, fqn: String, value: WdlValue) =>
      // Append this new output to the DB's event journal as a simple KVP (Output => {fqn, }
      // Also append to the specialised outputs table?

      // This replaces the functionality on WorkflowManagerActor - the slick endpoint would use this service actor instead.
      // NB: Move the messages WorkflowMetadata and WorkflowMetadataResponse into the MetadataServiceActor object instead.
    case WorkflowMetadata(workflowId) =>
      val workflowMetadataResponse: WorkflowMetadataResponse = ??? // Build it from the database
      sender ! workflowMetadataResponse
  }

}
