package wdl.draft3.transforms

import common.transforms.CheckedAtoB
import wdl.model.draft3.elements.ExpressionElement.IdentifierLookup
import wdl.model.draft3.elements.{DeclarationElement, ExpressionElement, InputDeclarationElement}
import wdl.draft3.transforms.wdlom2wom.ValueConsumer.ops._

package object wdlom2wom {
  val workflowDefinitionElementToWomWorkflowDefinition = CheckedAtoB.fromErrorOr(WorkflowDefinitionElementToWomWorkflowDefinition.convert)
  val fileElementToWomExecutable = CheckedAtoB.fromErrorOr(FileElementToWomExecutable.convert)

  implicit object expressionElementValueConsumer extends ValueConsumer[ExpressionElement] {
    override def consumedValues(a: ExpressionElement): Set[ConsumedValue] = a match {
      case IdentifierLookup(id) => Set(ConsumedSingleValue(id))
      case _ => Set.empty // TODO fill in other expression types
    }
  }

  implicit object unlinkedGraphNodeValueGeneratorAndConsumer extends ValueGenerator[UnlinkedGraphNode] with ValueConsumer[UnlinkedGraphNode] {
    override def consumedValues(a: UnlinkedGraphNode): Set[ConsumedValue] = a.fromElement match {
      case InputDeclarationElement(_, _, None) => Set.empty
      case DeclarationElement(_, _, Some(expr)) => expr.consumedValues
    }

    override def generatedValues(a: UnlinkedGraphNode): Set[String] = a.fromElement match {
      case DeclarationElement(_, name, _) => Set(name)
    }
  }
}
