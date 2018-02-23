package wdl.draft3.transforms.wdlom2wom.graph

import wdl.draft3.transforms.wdlom2wom.expression.linking._
import wdl.draft3.transforms.wdlom2wom.linking.UnlinkedValueConsumer.ops._
import wdl.draft3.transforms.wdlom2wom.linking._
import wdl.model.draft3.elements.{DeclarationElement, InputDeclarationElement, WorkflowGraphElement}

package object linking {
  implicit object graphElementUnlinkedValueGenerator extends UnlinkedValueGenerator[WorkflowGraphElement] {

    override def generatedValueNames(a: WorkflowGraphElement): Set[UnlinkedGeneratedValueName] = a match {
      case DeclarationElement(_, name, _) => Set(UnlinkedIdentifierName(name))
      // TODO fill in other expression types
      case other => throw new Exception(s"Cannot generate consumed values for WorkflowGraphNodeElement $other")
    }
  }

  implicit object graphElementUnlinkedValueConsumer extends UnlinkedValueConsumer[WorkflowGraphElement] {
    override def consumedValueNames(a: WorkflowGraphElement): Set[UnlinkedConsumedValueName] = a match {
      case InputDeclarationElement(_, _, None) => Set.empty
      case DeclarationElement(_, _, Some(expr)) => expr.consumedValueNames
      // TODO fill in other expression types
      case other => throw new Exception(s"Cannot generate consumed values for WorkflowGraphNodeElement $other")
    }
  }
}
