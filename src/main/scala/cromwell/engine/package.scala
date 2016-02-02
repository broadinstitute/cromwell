package cromwell

import java.nio.file.{Path, Paths}

import cromwell.engine.io.gcs.GcsFileSystem
import org.joda.time.DateTime
import org.slf4j.Logger
import wdl4s._
import wdl4s.values.{WdlFile, WdlValue}

import scala.language.implicitConversions
import scala.util.{Failure, Try}
import scalaz.ValidationNel

package object engine {

  class WorkflowContext(val root: String)
  class CallContext(override val root: String, val stdout: String, val stderr: String) extends WorkflowContext(root)

  /**
   * Represents the collection of source files that a user submits to run a workflow
   */
  final case class WorkflowSourceFiles(wdlSource: WdlSource, inputsJson: WdlJson, workflowOptionsJson: WorkflowOptionsJson)

  final case class AbortFunction(function: ()=>Unit)
  final case class AbortRegistrationFunction(register: AbortFunction=>Unit)

  final case class ExecutionEventEntry(description: String, startTime: DateTime, endTime: DateTime)
  final case class ExecutionHash(overallHash: String, dockerHash: Option[String])

  type ErrorOr[+A] = ValidationNel[String, A]

  case class SymbolHash(value: String) extends Ordered[SymbolHash] {
    def compare(that: SymbolHash) = this.value compare that.value
  }
  type FileHasher = WdlFile => SymbolHash

  type WorkflowOptionsJson = String
  type WorkflowOutputs = Map[FullyQualifiedName, CallOutput]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  case class CallOutput(wdlValue: WdlValue, hash: Option[SymbolHash])
  type CallOutputs = Map[LocallyQualifiedName, CallOutput]
  type HostInputs = Map[String, WdlValue]

  class CromwellFatalException(exception: Throwable) extends Exception(exception)

  implicit class EnhancedFullyQualifiedName(val fqn: FullyQualifiedName) extends AnyVal {
    def isScatter = fqn.contains(Scatter.FQNIdentifier)
    def scopeAndVariableName: (String, String) = {
      val array = fqn.split("\\.(?=[^\\.]+$)")
      (array(0), array(1))
    }
  }

  implicit class EnhancedCallOutputMap[A](val m: Map[A, CallOutput]) extends AnyVal {
    def mapToValues: Map[A, WdlValue] = m map {
      case (k, CallOutput(wdlValue, hash)) => (k, wdlValue)
    }
  }

  object PathString {
    implicit class UriString(val str: String) extends AnyVal {
      def isGcsUrl: Boolean = str.startsWith("gs://")

      def isUriWithProtocol: Boolean = "^[a-z]+://".r.findFirstIn(str).nonEmpty

      def toPath(workflowLogger: Logger, gcsFileSystem: Try[GcsFileSystem] = Failure(new Throwable("No GCS Filesystem"))): Path = {
        str match {
          case path if path.isGcsUrl && gcsFileSystem.isSuccess => gcsFileSystem.get.getPath(str)
          case path if path.isGcsUrl => throw new Throwable(s"Unable to parse GCS path $path: ${gcsFileSystem.failed.get.getMessage}")
          case path if !path.isUriWithProtocol => Paths.get(path)
          case path => throw new Throwable(s"Unable to parse $path")
        }
      }
    }
  }
}
