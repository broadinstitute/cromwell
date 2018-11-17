package cromwell.engine.instrumentation

import akka.actor.Actor
import cats.data.NonEmptyList
import cromwell.core.WorkflowState
import cromwell.core.instrumentation.InstrumentationPrefixes._
import cromwell.engine.instrumentation.WorkflowInstrumentation._
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowStoreState
import cromwell.services.instrumentation.CromwellInstrumentation._
import cromwell.services.instrumentation.CromwellInstrumentationActor

import scala.concurrent.duration.FiniteDuration
import scala.language.postfixOps

object WorkflowInstrumentation {
  private val WorkflowStatePaths: Map[WorkflowState, InstrumentationPath] = WorkflowState.WorkflowStateValues map { state =>
    state -> NonEmptyList.of(state.toString)
  } toMap

  // Use "Queued" instead of "Submitted" as it seems to reflect better the actual state
  private val SubmittedPath = NonEmptyList.of("Queued")
  private val RunningPath = NonEmptyList.of(WorkflowStoreState.Running.toString)
  private val OnHoldPath = NonEmptyList.of(WorkflowStoreState.OnHold.toString)
  private val AbortingPath = NonEmptyList.of(WorkflowStoreState.Aborting.toString)
}

/**
  * Provides helper methods for workflow instrumentation
  */
trait WorkflowInstrumentation extends CromwellInstrumentationActor { this: Actor =>
  private def workflowStatePath(workflowState: WorkflowState): InstrumentationPath = WorkflowInstrumentation.WorkflowStatePaths(workflowState)

  /**
    * Generic method to increment a workflow related counter metric value
    */
  def incrementWorkflow(statsDPath: InstrumentationPath): Unit = {
    increment(statsDPath, WorkflowPrefix)
  }

  /**
    * Generic method to add a workflow related timing metric value
    */
  def setTimingWorkflow(statsDPath: InstrumentationPath, duration: FiniteDuration): Unit = {
    sendTiming(statsDPath, duration, WorkflowPrefix)
  }

  /**
    * Generic method to update a workflow related gauge metric value
    */
  def sendGaugeWorkflow(statsDPath: InstrumentationPath, value: Long): Unit = {
    sendGauge(statsDPath, value, WorkflowPrefix)
  }

  /**
    * Counts every time a workflow enters a given state
    */
  def incrementWorkflowState(workflowState: WorkflowState): Unit = {
    incrementWorkflow(workflowStatePath(workflowState))
  }

  /**
    * Add a timing value for the run time of a workflow in a given state
    */
  //* TODO: enforce a terminal state ?
  def setWorkflowTimePerState(workflowState: WorkflowState, duration: FiniteDuration): Unit = {
    setTimingWorkflow(workflowStatePath(workflowState), duration)
  }

  /**
    * Set the current number of submitted workflows (queued but not running)
    */
  def updateWorkflowsQueued(count: Int) = sendGaugeWorkflow(SubmittedPath, count.toLong)

  /**
    * Set the current number of running workflows
    */
  def updateWorkflowsRunning(count: Int) = sendGaugeWorkflow(RunningPath, count.toLong)

  /**
    * Set the current number of on hold workflows
    */
  def updateWorkflowsOnHold(count: Int) = sendGaugeWorkflow(OnHoldPath, count.toLong)

  /**
    * Set the current number of aborting workflows
    */
  def updateWorkflowsAborting(count: Int) = sendGaugeWorkflow(AbortingPath, count.toLong)
}
