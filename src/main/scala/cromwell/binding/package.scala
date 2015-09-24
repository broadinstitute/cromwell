package cromwell

import com.typesafe.config.ConfigFactory
import spray.json._
import scala.util.{Try, Success, Failure}
import cromwell.binding.values.WdlValue
import cromwell.engine.WorkflowId
import cromwell.parser.BackendType

/**
 * ==WDL Bindings for Scala==
 *
 * This package contains an API to convert WDL ASTs into native Scala objects.
 *
 * Clients of this package should be able to interface with a WDL file without ever coming in direct contact with an Abstract Syntax Tree
 *
 * Internally, this package is built on top of [[cromwell.parser]].
 */

package object binding {
  type WdlSource = String
  type WdlJson = String
  type WorkflowOptionsJson = String
  type WorkflowOptions = Map[String, String]
  type WorkflowRawInputs = Map[FullyQualifiedName, Any]
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WdlValue]
  type WorkflowOutputs = Map[FullyQualifiedName, WdlValue]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type CallInputs = Map[String, WdlValue]
  type CallOutputs = Map[LocallyQualifiedName, WdlValue]
  type HostInputs = Map[String, WdlValue]
  type ImportResolver = String => WdlSource

  /**
   * Constructs a representation of a particular workflow invocation.  As with other
   * case classes and apply() methods, this will throw an exception if it cannot be
   * created
   */
  case class WorkflowDescriptor(id: WorkflowId, sourceFiles: WorkflowSourceFiles) {
    val workflowOptions = Try(sourceFiles.workflowOptionsJson.parseJson) match {
      case Success(JsObject(options)) =>
        if (options.values.nonEmpty && !options.values.exists(_.isInstanceOf[JsString])) {
          throw new Throwable(s"Workflow ${id.toString} options JSON is not a String -> String map: ${sourceFiles.workflowOptionsJson}")
        }
        options map { case (k, v) => k -> v.asInstanceOf[JsString].value }
      case _ => throw new Throwable(s"Workflow ${id.toString} contains bad workflow options JSON: ${sourceFiles.inputsJson}")
    }

    val backendType = BackendType.from(workflowOptions.getOrElse("default_backend", ConfigFactory.load.getConfig("backend").getString("backend")))
    val namespace = NamespaceWithWorkflow.load(sourceFiles.wdlSource, backendType)
    val name = namespace.workflow.name
    val shortId = id.toString.split("-")(0)

    val rawInputs = Try(sourceFiles.inputsJson.parseJson) match {
      case Success(JsObject(inputs)) => inputs
      case _ => throw new Throwable(s"Workflow ${id.toString} contains bad inputs JSON: ${sourceFiles.inputsJson}")
    }

    // Currently we are throwing an exception if construction of the workflow descriptor fails, hence .get on the Trys
    val coercedInputs = namespace.coerceRawInputs(rawInputs).get
    val declarations = namespace.staticDeclarationsRecursive(coercedInputs).get
    val actualInputs: WorkflowCoercedInputs = coercedInputs ++ declarations
  }

  case class WorkflowSourceFiles(wdlSource: WdlSource, inputsJson: WdlJson, workflowOptionsJson: WorkflowOptionsJson)

  implicit class EnhancedFullyQualifiedName(val fqn: FullyQualifiedName) extends AnyVal {
    def isScatter = fqn.contains(Scatter.FQNIdentifier)
  }
}
