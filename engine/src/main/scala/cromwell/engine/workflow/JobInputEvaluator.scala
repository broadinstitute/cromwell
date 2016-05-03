package cromwell.engine.workflow

import cromwell.backend.BackendJobDescriptorKey
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.ExecutionIndex._
import cromwell.engine.workflow.lifecycle.WorkflowExecutionActor.{OutputCallKey, OutputEntry, OutputStore}
import cromwell.util.TryUtil
import wdl4s._
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values.{WdlCallOutputsObject, WdlValue}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

/**
   * Resolve and evaluate call input declarations
   *
   * As a result, the following wdl
   *
   * task t {
   *   File q
   *   String r = read_string(q) #(1)
   *   String s
   * }
   *
   * workflow w {
   *   File f
   *   String a = read_string(f) #(2)
   *   call t { input:
   *              s = a + read_string(f), #(3)
   *              q = f
   *          }
   * }
   *
   * will be evaluated such that
   *
   * read_string (1) will be evaluated using BFN
   * read_string (2) will be evaluated using EFN
   * read_string (3) will be evaluated using BFN
   *
   * Where
   *   BFN is the WDL function implementation from the backend t will run on
   *   EFN is the WDL function implementation from the engine
   *
   * q, r and s will be evaluated sequentially, such that when evaluating r,
   * the value of q has already been resolved in the appropriate scope and evaluated by the appropriate wdl functions.
   */
class JobInputEvaluator(workflowDescriptor: EngineWorkflowDescriptor) {

  private val calls = workflowDescriptor.backendDescriptor.workflowNamespace.workflow.calls

  /*
   *  /|\ this could be very time consuming, as it potentially involves IO when using engine functions
   * /_._\
   */
  def resolveAndEvaluate(jobKey: BackendJobDescriptorKey,
                         wdlFunctions: WdlStandardLibraryFunctions,
                         outputStore: OutputStore): Try[Map[LocallyQualifiedName, WdlValue]] = {
    val call = jobKey.call
    lazy val callInputsFromFile = unqualifiedInputsFromInputFile(call)
    lazy val workflowScopedLookup = buildWorkflowScopedLookup(jobKey, outputStore) _

    // Try to resolve, evaluate and coerce declarations in order
    val inputEvaluationAttempt = call.task.declarations.foldLeft(Map.empty[LocallyQualifiedName, Try[WdlValue]])((inputs, declaration) => {
      val name = declaration.name

      // Try to resolve the declaration, and upon success evaluate the expression
      // If the declaration is resolved but can't be evaluated this will throw an evaluation exception
      // If it can't be resolved it's ignored and won't appear in the final input map
      val evaluated: Option[Try[WdlValue]] = declaration.expression match {
        // Static expression in the declaration
        case Some(expr) => Option(expr.evaluate(buildMapBasedLookup(inputs), wdlFunctions))
        // Expression found in the input mappings
        case None if call.inputMappings.contains(name) => Option(call.inputMappings(name).evaluate(workflowScopedLookup, wdlFunctions))
        // Expression found in the input file
        case None if callInputsFromFile.contains(name) => Option(Success(callInputsFromFile(name)))
        // Expression can't be found
        case _ => None
      }

      // Leave out unresolved declarations
      evaluated match {
        case Some(value) =>
          val coercedValue = value flatMap declaration.wdlType.coerceRawValue
          inputs + ((name, coercedValue))
        case None => inputs
      }
    })

    TryUtil.sequenceMap(inputEvaluationAttempt, s"Input evaluation for Call ${call.fullyQualifiedName} failed")
  }

  private def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
    val lastIndex = fullyQualifiedName.lastIndexOf(".")
    (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
  }

  // Split the fqn value into 2 chunks, separating the last member out. eg: "myWorkflow.myTask.myInput" => ("myWorkflow.myTask", "myInput")
  private lazy val splitInputs = workflowDescriptor.backendDescriptor.inputs map {
    case (fqn, v) => splitFqn(fqn) -> v
  }

  // Unqualified workflow level inputs
  private lazy val unqualifiedWorkflowInputs: Map[LocallyQualifiedName, WdlValue] = splitInputs collect {
    case((root, inputName), v) if root == workflowDescriptor.namespace.workflow.unqualifiedName => inputName -> v
  }

  // Unqualified call inputs for a specific call, from the input json
  private def unqualifiedInputsFromInputFile(call: Call): Map[LocallyQualifiedName, WdlValue] = splitInputs collect {
    case((root, inputName), v) if root == call.fullyQualifiedName => inputName -> v
  }

  // In this case, the lookup function is effectively equivalent to looking into unqualifiedWorkflowInputs for the value
  // because the resolution / evaluation / coercion has already happened in the MaterializeWorkflowDescriptorActor
  private def lookupWorkflowDeclaration(identifier: String) = {
    unqualifiedWorkflowInputs.get(identifier) match {
      case Some(value) => Success(value)
      case None => Failure(new WdlExpressionException(s"Could not resolve variable $identifier as a workflow input"))
    }
  }

  private def lookupCall(key: BackendJobDescriptorKey, outputStore: OutputStore)(identifier: String): Try[WdlCallOutputsObject] = {
    calls find { _.unqualifiedName == identifier } match {
      case Some(matchedCall) =>
        /**
          * After matching the Call, this determines if the `key` depends on a single shard
          * of a scatter'd job or if it depends on the whole thing.  Right now, the heuristic
          * is "If we're both in a scatter block together, then I depend on a shard.  If not,
          * I depend on the collected value"
          *
          * TODO: nested-scatter - this will likely not be sufficient for nested scatters
          */
        val index: ExecutionIndex = matchedCall.closestCommonAncestor(key.scope) flatMap {
          case s: Scatter => key.index
          case _ => None
        }
        fetchCallOutputEntries(matchedCall, index, outputStore)
      case None => Failure(new WdlExpressionException(s"Could not find a call with identifier '$identifier'"))
    }
  }

  private def fetchCallOutputEntries(call: Call, index: ExecutionIndex, outputStore: OutputStore): Try[WdlCallOutputsObject] = {
    def outputEntriesToMap(outputs: Traversable[OutputEntry]): Map[String, Try[WdlValue]] = {
      outputs map { output =>
        output.wdlValue match {
          case Some(wdlValue) => output.name -> Success(wdlValue)
          case None => output.name -> Failure(new RuntimeException(s"Could not retrieve output ${output.name} value"))
        }
      } toMap
    }

    outputStore.get(OutputCallKey(call, index)) match {
      case Some(outputs) =>
        TryUtil.sequenceMap(outputEntriesToMap(outputs), s"Output fetching for call ${call.unqualifiedName}") map { outputsMap =>
          WdlCallOutputsObject(call, outputsMap)
        }
      case None => Failure(new RuntimeException(s"Could not find call ${call.unqualifiedName}"))
    }
  }

  private def buildWorkflowScopedLookup(jobKey: BackendJobDescriptorKey, outputStore: OutputStore)(identifier: String): WdlValue = {
    // Will need to be enhanced with scatter variable resolvers
    val resolvers = Stream(lookupWorkflowDeclaration _, lookupCall(jobKey, outputStore) _)

    resolvers map { _(identifier) } find { _.isSuccess } getOrElse {
      Failure(new WdlExpressionException(s"Could not resolve variable $identifier"))
    } get
  }

  private def buildMapBasedLookup(evaluatedDeclarations: Map[LocallyQualifiedName, Try[WdlValue]])(identifier: String): WdlValue = {
    val successfulEvaluations = evaluatedDeclarations collect {
      case (k, v) if v.isSuccess => k -> v.get
    }
    successfulEvaluations.getOrElse(identifier, throw new WdlExpressionException(s"Could not resolve variable $identifier as a task input"))
  }
}
