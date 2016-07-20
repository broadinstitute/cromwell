package cromwell.engine.workflow

import cromwell.core.{WorkflowId, WorkflowSourceFiles}

import scala.concurrent.{ExecutionContext, Future}
import scalaz.NonEmptyList

package object workflowstore {

  sealed trait WorkflowState {def isStartable: Boolean}

  case object Running extends WorkflowState { override def isStartable = false }
  sealed trait StartableState extends WorkflowState { override def isStartable = true }
  case object Submitted extends StartableState
  case object Restartable extends StartableState

  final case class WorkflowToStart(id: WorkflowId, sources: WorkflowSourceFiles, state: StartableState)
}
