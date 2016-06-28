package cromwell.backend.impl.jes

import akka.actor.{ActorRef, Props}
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend._
import cromwell.backend.async.AsyncBackendJobExecutionActor.Execute
import org.slf4j.LoggerFactory

import scala.concurrent.{Future, Promise}
import scala.language.postfixOps

object JesJobExecutionActor {
  val logger = LoggerFactory.getLogger("JesBackend")

  def props(jobDescriptor: BackendJobDescriptor, jesWorkflowInfo: JesConfiguration, initializationData: JesBackendInitializationData): Props = {
    Props(new JesJobExecutionActor(jobDescriptor, jesWorkflowInfo, initializationData))
  }
}

case class JesJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                jesConfiguration: JesConfiguration,
                                initializationData: JesBackendInitializationData)
  extends BackendJobExecutionActor {

  override def receive: Receive = LoggingReceive {
    case AbortJobCommand =>
      executor.foreach(_ ! AbortJobCommand)
    case abortResponse: AbortedResponse =>
      context.parent ! abortResponse
      context.stop(self)

    // PBE TODO: use PartialFunction.orElse instead of this catch-all case, because Akka might use isDefinedAt over this partial function
    case message => super.receive(message)
  }

  override val configurationDescriptor = jesConfiguration.configurationDescriptor

  // PBE keep a reference to be able to hand to a successor executor if the failure is recoverable, or to complete as a
  // failure if the failure is not recoverable.
  private lazy val completionPromise = Promise[BackendJobExecutionResponse]()

  // PBE keep a reference for abort purposes.  Maybe we don't need this if we can look up actors by name.
  private var executor: Option[ActorRef] = None

  // PBE there should be some consideration of supervision here.

  override def recover: Future[BackendJobExecutionResponse] = ???

  override def execute: Future[BackendJobExecutionResponse] = {
    val executorRef = context.actorOf(JesAsyncBackendJobExecutionActor.props(jobDescriptor, completionPromise, jesConfiguration, initializationData))
    executor = Option(executorRef)
    executorRef ! Execute
    completionPromise.future
  }

  override def abort: Unit = {}
}
