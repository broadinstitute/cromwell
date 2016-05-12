package cromwell.engine.workflow.profiler

import akka.actor._
import akka.event.LoggingReceive
import com.typesafe.scalalogging.StrictLogging
import cromwell.core.WorkflowId
import cromwell.engine.WorkflowState
import cromwell.engine.workflow.WorkflowActor.{InitializingWorkflowState, WorkflowActorStateAndData}
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.{WorkflowExecutionActorStateAndData, WorkflowExecutionSuccessfulState}
import cromwell.engine.workflow.{CromwellFsmStateAndData, WorkflowMetadataKeys}
import cromwell.services.MetadataServiceActor._
import org.joda.time.DateTime

import scala.concurrent.{ExecutionContext, Future}

object WorkflowProfilerActor {

  case class WorkflowManagerMetadataMessage(workflowId: String, submissionTime: String)

  implicit class EnhancedServiceRegistry(val serviceRegistryActor: ActorRef) extends StrictLogging {
    def publishToMetadataService(message: MetadataEvent)(implicit actorSystem: ActorSystem, ec: ExecutionContext): Future[MetadataServiceResponse] = {
      import cromwell.util.PromiseActor._
      // Since this will always publish, we expect the callee to just provide an event,
      // and we wrap it into a `PutMetadataAction` before sending it to the MetadataService
      serviceRegistryActor.askNoTimeout(PutMetadataAction(message)).mapTo[MetadataServiceResponse].recover {
        case reason: Throwable =>
          // TODO: Do we still proceed with the execution? What do we do with the exception here?
          logger.warn(s"Failed to push event $message in the MetadataService: ${reason.getMessage}", reason)
          throw reason
      }
    }
  }

  def props(workflowId: WorkflowId, serviceRegistryActor: ActorRef): Props = {
    Props(new WorkflowProfilerActor(workflowId, serviceRegistryActor))
  }
}

case class WorkflowProfilerActor(workflowId: WorkflowId, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {

  import WorkflowProfilerActor._

  private val tag = s"WorkflowProfiler-${workflowId.toString}"
  implicit val actorSystem = context.system
  implicit val ec = context.dispatcher

  override def preStart(): Unit = {
    // Subscribe to the events published on the stream
    context.system.eventStream.subscribe(self, classOf[CromwellFsmStateAndData])
  }

  override def receive = LoggingReceive {
    case workflowManagerActorMessage: WorkflowManagerMetadataMessage => handleMessageFromWorkflowManagerActor(workflowManagerActorMessage)
    case workflowActorMessage: WorkflowActorStateAndData => handleMessageFromWorkflowActor(workflowActorMessage)
    case workflowExecutionActorMessage: WorkflowExecutionActorStateAndData => handleMessageFromWorkflowExecutionActor(workflowExecutionActorMessage)
    case unhandledMessage: Any => log.warning(s"$tag Unknown Message received: $unhandledMessage")
  }

  private def handleMessageFromWorkflowManagerActor(workflowManagerMetadataMessage: WorkflowManagerMetadataMessage): Unit = {
    val metadataEventMsgs = List(
      MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Id), MetadataValue(workflowManagerMetadataMessage.workflowId)),
      MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.SubmissionTime), MetadataValue(workflowManagerMetadataMessage.submissionTime)),
      // Currently, submission time is the same as start time
      MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.StartTime), MetadataValue(workflowManagerMetadataMessage.submissionTime))
    )
    metadataEventMsgs foreach serviceRegistryActor.publishToMetadataService
  }

  private def handleMessageFromWorkflowActor(workflowActorStateAndData: WorkflowActorStateAndData): Unit = {
    workflowActorStateAndData match {
      case WorkflowActorStateAndData(state, data) if state == InitializingWorkflowState =>
        // `get` is safe since in this state we are guaranteed to have a ValidWorkflowDescriptor
        val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Name), MetadataValue(data.workflowDescriptor.get.name))
        serviceRegistryActor.publishToMetadataService(metadataEventMsg)
      case WorkflowActorStateAndData(state, data) =>
        // This updates the workflow status
        publishCurrentStateToMetadataService(state.workflowState)
        if (state.terminal) {
          // Add the end time of the workflow in the MetadataService
          val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.EndTime), MetadataValue(DateTime.now.toString))
          serviceRegistryActor.publishToMetadataService(metadataEventMsg)
          // As the state is terminal, this is expected to be the last message related to a particular Workflow. This actor's job here seems to be done
          self ! PoisonPill
        }
    }
  }

  private def handleMessageFromWorkflowExecutionActor(workflowExecutionActorStateAndData: WorkflowExecutionActorStateAndData) = {
    workflowExecutionActorStateAndData match {
      case WorkflowExecutionActorStateAndData(state, data) if state == WorkflowExecutionSuccessfulState =>
        val outputs = data.outputStore.store.values.flatten
        val keyValues = data.outputStore.store.flatMap{
          case (key, value) => value.map(entry => s"${key.call.fullyQualifiedName}.${entry.name}" -> entry.wdlValue.map(_.toWdlString).getOrElse("NA"))
        }
        val metadataEventMsgs = keyValues.map{ case (k,v) => MetadataEvent(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Outputs}.$k"), MetadataValue(v))}
        metadataEventMsgs foreach serviceRegistryActor.publishToMetadataService
    }
  }

  /* Update the current State of the Workflow (corresponding to the FSM state) in the Metadata service */
  private def publishCurrentStateToMetadataService(workflowState: WorkflowState): Unit = {
    val metadataEventMsg = MetadataEvent(MetadataKey(workflowId, None, WorkflowMetadataKeys.Status), MetadataValue(workflowState.toString))
    serviceRegistryActor.publishToMetadataService(metadataEventMsg)
  }
}
