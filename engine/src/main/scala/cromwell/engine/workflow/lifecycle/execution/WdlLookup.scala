package cromwell.engine.workflow.lifecycle.execution

import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.ExecutionIndex._
import wdl4s._
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values.{WdlArray, WdlCallOutputsObject, WdlValue}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

trait WdlLookup {

  def workflowDescriptor: EngineWorkflowDescriptor
  def executionStore: ExecutionStore
  def outputStore: OutputStore
  def expressionLanguageFunctions: WdlStandardLibraryFunctions

  private lazy val splitInputs = workflowDescriptor.backendDescriptor.inputs map {
    case (fqn, v) => splitFqn(fqn) -> v
  }

  // Unqualified workflow level inputs
  private lazy val unqualifiedWorkflowInputs: Map[LocallyQualifiedName, WdlValue] = splitInputs collect {
    case((root, inputName), v) if root == workflowDescriptor.namespace.workflow.unqualifiedName => inputName -> v
  }

  /**
    * Lookup an identifier by
    * first looking at the completed calls map
    * and if not found traversing up the scope hierarchy from the scope from which the lookup originated.
    */
  def hierarchicalLookup(scope: Scope, index: ExecutionIndex)(identifier: String): WdlValue = {
    // First lookup calls
    lookupCall(scope, index, identifier) recoverWith {
      // Lookup in the same scope (currently no scope support this but say we have scatter declarations, or multiple scatter variables, or nested workflows..)
      case _: VariableNotFoundException | _: WdlExpressionException => scopedLookup(scope, index, identifier)
    } recover {
      // Lookup parent if present
      case _: VariableNotFoundException | _: WdlExpressionException => scope.parent match {
        case Some(parent) => hierarchicalLookup(parent, index)(identifier)
        case None => throw new VariableNotFoundException(s"Can't find $identifier")
      }
    } get
  }

  private def scopedLookup(scope: Scope, index: ExecutionIndex, identifier: String): Try[WdlValue] = {
    def scopedLookupFunction = scope match {
      case scatter: Scatter if index.isDefined => lookupScatter(scatter, index.get) _
      case workflow: Workflow => lookupWorkflowDeclaration _
      case _ => (_: String) => Failure(new VariableNotFoundException(s"Can't find $identifier in scope $scope"))
    }

    scopedLookupFunction(identifier)
  }

  // In this case, the scopedLookup function is effectively equivalent to looking into unqualifiedWorkflowInputs for the value
  // because the resolution / evaluation / coercion has already happened in the MaterializeWorkflowDescriptorActor
  private def lookupWorkflowDeclaration(identifier: String) = {
    unqualifiedWorkflowInputs.get(identifier) match {
      case Some(value) => Success(value)
      case None => Failure(new WdlExpressionException(s"Could not resolve variable $identifier as a workflow input"))
    }
  }

  private def lookupScatter(scatter: Scatter, index: Int)(identifier: String): Try[WdlValue] = {
    if (identifier == scatter.item) {
      // Scatters are not indexed yet (they can't be nested)
      val scatterLookup = hierarchicalLookup(scatter, None) _
      scatter.collection.evaluate(scatterLookup, expressionLanguageFunctions) map {
        case collection: WdlArray if collection.value.isDefinedAt(index) => collection.value(index)
        case collection: WdlArray => throw new RuntimeException(s"Index $index out of bound in $collection for scatter ${scatter.fullyQualifiedName}")
        case other => throw new RuntimeException(s"Scatter ${scatter.fullyQualifiedName} collection is not an array: $other")
      } recover {
        case e => throw new RuntimeException(s"Failed to evaluate collection for scatter ${scatter.fullyQualifiedName}", e)
      }
    } else {
      Failure(new VariableNotFoundException(identifier))
    }
  }

  private def lookupCall(scope: Scope, scopeIndex: ExecutionIndex, identifier: String): Try[WdlCallOutputsObject] = {
    val calls = executionStore.store.keys.view map { _.scope } collect { case c: Call => c }

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
        val index: ExecutionIndex = matchedCall.closestCommonAncestor(scope) flatMap {
          case s: Scatter => scopeIndex
          case _ => None
        }

        outputStore.fetchCallOutputEntries(matchedCall, index)
      case None => Failure(new WdlExpressionException(s"Could not find a call with identifier '$identifier'"))
    }
  }

}
