package wdl.draft3.transforms.linking

import cats.syntax.validated._

import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.graph.UnlinkedValueConsumer.ops._
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.draft3.transforms.linking.expression._
import wdl.draft3.transforms.linking.typemakers._
import wdl.model.draft3.elements.{DeclarationElement, InputDeclarationElement, WorkflowGraphElement}
import wdl.model.draft3.graph._
import wom.types.WomType

package object graph {
  implicit object graphElementUnlinkedValueGenerator extends UnlinkedValueGenerator[WorkflowGraphElement] {

    override def generatedValueHandles(a: WorkflowGraphElement, typeAliases: Map[String, WomType]): ErrorOr[Set[GeneratedValueHandle]] = a match {
      case DeclarationElement(typeElement, name, _) =>
        typeElement.determineWomType(typeAliases) map { t => Set(GeneratedIdentifierValueHandle(name, t)) }

      // TODO fill in other expression types
      case other => s"Cannot generate consumed values for WorkflowGraphNodeElement $other".invalidNel
    }
  }

  implicit object graphElementUnlinkedValueConsumer extends UnlinkedValueConsumer[WorkflowGraphElement] {
    override def consumedValueHooks(a: WorkflowGraphElement): Set[UnlinkedConsumedValueHook] = a match {
      case InputDeclarationElement(_, _, None) => Set.empty
      case DeclarationElement(_, _, Some(expr)) => expr.consumedValueHooks
      // TODO fill in other expression types
      case other => throw new Exception(s"Cannot generate consumed values for WorkflowGraphNodeElement $other")
    }
  }
}
