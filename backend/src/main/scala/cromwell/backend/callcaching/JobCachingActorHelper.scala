package cromwell.backend.callcaching

import akka.actor.Actor
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import cromwell.core.logging.JobLogging

/**
  * Base trait for mixins related to caching functionality.
  *
  * This trait is extended by sub traits that provide more functionality for a specific backend. That sub trait  is
  * then mixed into a pair of actors, the first that does the (async-)execution, the second that responds to cache hits.
  *
  * Very closely tied to, but not the same as a BackendLifecycleActor/BackendJobLifecycleActor. Those traits define and
  * implement other methods. Plus, the async-executor actors do not actually extend from BackendJobLifecycleActor.
  *
  * Child job caching actor helper trait implementations should also include the self type:
  * {{{
  *   this: Actor with JobLogging =>
  * }}}
  */
trait JobCachingActorHelper {
  this: Actor with JobLogging =>
  // For Logging and boilerplate
  override lazy val workflowId = jobDescriptor.workflowDescriptor.id
  override lazy val jobTag = jobDescriptor.key.tag

  /**
    * The job being executed or cache hit.
    */
  def jobDescriptor: BackendJobDescriptor

  /**
    * The configuration for the backend, in the context of the entire Cromwell configuration file.
    */
  def configurationDescriptor: BackendConfigurationDescriptor
}
