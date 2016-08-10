package cromwell.backend.impl.jes

import akka.actor.{ActorRef, Props}
import akka.event.LoggingReceive
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend._
import cromwell.backend.async.AsyncBackendJobExecutionActor.{Execute, Recover}
import cromwell.backend.impl.jes.JesAsyncBackendJobExecutionActor.JesJobId
import org.slf4j.LoggerFactory

import scala.concurrent.{Future, Promise}
import scala.language.postfixOps
import cromwell.backend.impl.jes.JesJobExecutionActor._
import cromwell.services.keyvalue.KeyValueServiceActor._

object JesJobExecutionActor {
  val logger = LoggerFactory.getLogger("JesBackend")

  def props(jobDescriptor: BackendJobDescriptor,
            jesWorkflowInfo: JesConfiguration,
            initializationData: JesBackendInitializationData,
            serviceRegistryActor: ActorRef): Props = {
    Props(new JesJobExecutionActor(jobDescriptor, jesWorkflowInfo, initializationData, serviceRegistryActor))
  }

  val JesOperationIdKey = "__jes_operation_id"
}

case class JesJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                jesConfiguration: JesConfiguration,
                                initializationData: JesBackendInitializationData,
                                serviceRegistryActor: ActorRef)
  extends BackendJobExecutionActor {

  private def jesReceiveBehavior: Receive = LoggingReceive {
    case AbortJobCommand =>
      executor.foreach(_ ! AbortJobCommand)
    case abortResponse: AbortedResponse =>
      context.parent ! abortResponse
      context.stop(self)
    case KvPair(key, id @ Some(operationId)) if key.key == JesOperationIdKey =>
      // Successful operation ID lookup during recover.
      executor foreach { _ ! Recover(JesJobId(operationId))}
    case KvKeyLookupFailed(_) =>
      // Missed operation ID lookup during recover, fall back to execute.
      executor foreach { _ ! Execute }
    case KvFailure(_, e) =>
      // Failed operation ID lookup during recover, crash and let the supervisor deal with it.
      completionPromise.tryFailure(e)
      throw new RuntimeException("Failure attempting to look up JES operation ID for key " + jobDescriptor.key, e)
  }

  override def receive = jesReceiveBehavior orElse super.receive

  override val configurationDescriptor = jesConfiguration.configurationDescriptor

  private lazy val completionPromise = Promise[BackendJobExecutionResponse]()

  private var executor: Option[ActorRef] = None

  private def launchExecutor: Future[Unit] = Future {
    val executionProps = JesAsyncBackendJobExecutionActor.props(jobDescriptor,
      completionPromise,
      jesConfiguration,
      initializationData,
      serviceRegistryActor)
    val executorRef = context.actorOf(executionProps)
    executor = Option(executorRef)
    ()
  }

  override def recover: Future[BackendJobExecutionResponse] = {
    import JesJobExecutionActor._

    for {
      _ <- launchExecutor
      _ = serviceRegistryActor ! KvGet(ScopedKey(jobDescriptor.workflowDescriptor.id,
        KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt),
        JesOperationIdKey))
      c <- completionPromise.future
    } yield c
  }

  override def execute: Future[BackendJobExecutionResponse] = {
    for {
      _ <- launchExecutor
      _ = executor foreach { _ ! Execute }
      c <- completionPromise.future
    } yield c
  }

  override def abort: Unit = {}
}
