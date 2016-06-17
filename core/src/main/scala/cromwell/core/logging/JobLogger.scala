package cromwell.core.logging

import akka.actor.{Actor, ActorLogging}
import akka.event.LoggingAdapter
import cromwell.core.WorkflowId
import org.slf4j.Logger

trait JobLogging extends ActorLogging { this: Actor =>
  def workflowId: WorkflowId
  def jobTag: String

  lazy val jobLogger: Logger = new JobLogger(self.path.name, workflowId, jobTag, Option(log))
}

/**
  * Similar to WorkflowLogger, with the addition that the log tag includes the call tag.
  */
class JobLogger(name: String,
                workflowId: WorkflowId,
                jobTag: String,
                akkaLogger: Option[LoggingAdapter] = None,
                otherLoggers: Set[Logger] = Set.empty[Logger])
  extends WorkflowLogger(name, workflowId, akkaLogger, otherLoggers) {

  override def tag = s"$name [UUID(${workflowId.shortString})$jobTag]"
}
