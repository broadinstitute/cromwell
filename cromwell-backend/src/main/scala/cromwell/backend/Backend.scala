package cromwell.backend

import cromwell.backend.model.Subscription
import cromwell.caching.Cacheable

/**
  * Defines basic functionality to interact with a Backend.
  */
trait Backend extends Cacheable{
  /**
    * Prepare the task and context for execution.
    */
  def prepare(): Unit
  /**
    * Executes task in given context.
    */
  def execute(): Unit

  /**
    * Stops a task execution.
    */
  def stop(): Unit

  /**
    * Performs a cleanUp after the task was executed.
    */
  def cleanUp(): Unit

  /**
    * Subscribe to events on backend.
    */
  def subscribeToEvent[T](subscription: Subscription[T]): Unit

  /**
    * Unsubscribe to events on backend.
    */
  def unsubscribeToEvent[T](subscription: Subscription[T]): Unit
}