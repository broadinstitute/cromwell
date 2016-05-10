package cromwell.backend.impl.jes

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend._
import cromwell.backend.async.AsyncBackendJobExecutionActor.Execute
import org.slf4j.LoggerFactory

import scala.concurrent.{Future, Promise}
import scala.language.postfixOps


object JesJobExecutionActor {
  val logger = LoggerFactory.getLogger("JesBackend")

  def props(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor): Props =
    Props(new JesJobExecutionActor(jobDescriptor, configurationDescriptor))

}

case class JesJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                override val configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendJobExecutionActor {

  // PBE keep a reference to be able to hand to a successor executor if the failure is recoverable, or to complete as a
  // failure if the failure is not recoverable.
  private lazy val completionPromise = Promise[BackendJobExecutionResponse]()

  // PBE keep a reference for abort purposes.  Maybe we don't need this if we can look up actors by name.
  private var executor: Option[ActorRef] = None

  // PBE there should be some consideration of supervision here.

  /**
    * Restart or resume a previously-started job.
    */
  override def recover: Future[BackendJobExecutionResponse] = ???

  /**
    * Execute a new job.
    */
  override def execute: Future[BackendJobExecutionResponse] = {
    val executorRef = context.actorOf(JesAsyncBackendJobExecutionActor.props(jobDescriptor, configurationDescriptor, completionPromise))
    executor = Option(executorRef)
    executorRef ! Execute
    completionPromise.future
  }

  /**
    * Abort a running job.
    */
  override def abortJob: Unit = ???
}
