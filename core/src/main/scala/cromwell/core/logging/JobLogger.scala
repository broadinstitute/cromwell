package cromwell.core.logging

import akka.actor.{Actor, ActorLogging}
import akka.event.LoggingAdapter
import cromwell.core.{PossiblyNotRootWorkflowId, RootWorkflowId}
import org.slf4j.Logger

trait JobLogging extends ActorLogging { this: Actor =>
  def workflowIdForLogging: PossiblyNotRootWorkflowId
  def rootWorkflowIdForLogging: RootWorkflowId
  def jobTag: String

  lazy val jobLogger =
    new JobLogger(
      loggerName = self.path.name,
      workflowIdForLogging = workflowIdForLogging,
      rootWorkflowIdForLogging = rootWorkflowIdForLogging,
      jobTag = jobTag,
      akkaLogger = Option(log)
    )
}

/**
  * Similar to WorkflowLogger, with the addition that the log tag includes the call tag.
  */
class JobLogger(loggerName: String,
                workflowIdForLogging: PossiblyNotRootWorkflowId,
                rootWorkflowIdForLogging: RootWorkflowId,
                jobTag: String,
                akkaLogger: Option[LoggingAdapter] = None,
                otherLoggers: Set[Logger] = Set.empty[Logger])
  extends WorkflowLogger(
    loggerName = loggerName,
    workflowId = workflowIdForLogging,
    rootWorkflowId = rootWorkflowIdForLogging,
    akkaLogger = akkaLogger,
    otherLoggers = otherLoggers
  ) {

  override def tag = s"$loggerName [UUID(${workflowIdForLogging.shortString})$jobTag]"
}
