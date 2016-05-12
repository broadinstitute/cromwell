package cromwell

import com.typesafe.config.Config
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.joda.time.DateTime
import wdl4s._
import wdl4s.values.WdlValue

import scala.language.postfixOps
import scala.util.{Success, Try}

package object backend {

  trait JobKey {
    def scope: Scope
    def index: Option[Int]
    def attempt: Int
    def tag: String
  }

  /**
    * For uniquely identifying a job which has been or will be sent to the backend.
    */
  case class BackendJobDescriptorKey(call: Call, index: Option[Int], attempt: Int) extends JobKey {
    val scope = call
    private val indexString = index map { _.toString } getOrElse "NA"
    val tag = s"${call.fullyQualifiedName}:$indexString:$attempt"
    val isShard = index.isDefined
  }

  /**
    * For passing to a BackendWorkflowActor for job execution or recovery
    */
  case class BackendJobDescriptor(descriptor: BackendWorkflowDescriptor,
                                  key: BackendJobDescriptorKey,
                                  inputs: Map[LocallyQualifiedName, WdlValue]) {
    val call = key.call
  }

  /**
    * For passing to a BackendActor construction time
    */
  case class BackendWorkflowDescriptor(id: WorkflowId,
                                       workflowNamespace: NamespaceWithWorkflow,
                                       inputs: Map[FullyQualifiedName, WdlValue],
                                       workflowOptions: WorkflowOptions)

  /**
    * For passing to a BackendActor construction time
    */
  case class BackendConfigurationDescriptor(backendConfig: Config, globalConfig: Config)

  final case class ExecutionEventEntry(description: String, startTime: DateTime, endTime: DateTime)

  final case class ExecutionHash(overallHash: String, dockerHash: Option[String])

  final case class AttemptedLookupResult(name: String, value: Try[WdlValue]) {
    def toPair = name -> value
  }

  implicit class AugmentedAttemptedLookupSequence(s: Seq[AttemptedLookupResult]) {
    def toLookupMap: Map[String, WdlValue] = s collect {
      case AttemptedLookupResult(name, Success(value)) => (name, value)
    } toMap
  }

  case class PreemptedException(msg: String) extends Exception(msg)
}
