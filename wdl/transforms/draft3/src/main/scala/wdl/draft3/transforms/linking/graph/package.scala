package wdl.draft3.transforms.linking

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.graph.GraphElementValueConsumer
import wdl.model.draft3.graph.GraphElementValueConsumer.ops._
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.UnlinkedValueGenerator.ops._
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.draft3.transforms.linking.expression.consumed._
import wdl.draft3.transforms.linking.typemakers._
import wdl.model.draft3.elements.{DeclarationElement, _}
import wdl.model.draft3.graph._
import wom.callable.Callable
import wom.callable.Callable.OutputDefinition
import wom.types.{WomArrayType, WomOptionalType, WomType}

package object graph {
  implicit val graphElementUnlinkedValueGenerator: UnlinkedValueGenerator[WorkflowGraphElement] = new UnlinkedValueGenerator[WorkflowGraphElement] {

    override def generatedValueHandles(a: WorkflowGraphElement,
                                       typeAliases: Map[String, WomType],
                                       callables: Map[String, Callable]): ErrorOr[Set[GeneratedValueHandle]] = a match {
      case DeclarationElement(typeElement, name, _) =>
        typeElement.determineWomType(typeAliases) map { t => Set(GeneratedIdentifierValueHandle(name, t)) }
      case a: ScatterElement => a.generatedValueHandles(typeAliases, callables)
      case a: IfElement => a.generatedValueHandles(typeAliases, callables)
      case a: CallElement => a.generatedValueHandles(typeAliases, callables)
      case other => s"Cannot generate generated values for WorkflowGraphNodeElement $other".invalidNel
    }
  }

  implicit val scatterElementUnlinkedValueGenerator: UnlinkedValueGenerator[ScatterElement] = new UnlinkedValueGenerator[ScatterElement] {
    override def generatedValueHandles(a: ScatterElement, typeAliases: Map[String, WomType], callables: Map[String, Callable]): ErrorOr[Set[GeneratedValueHandle]] = {
      a.graphElements.toList.traverse(_.generatedValueHandles(typeAliases, callables)).map(_.toSet.flatten) map { _.map {
        case GeneratedIdentifierValueHandle(id, womType) => GeneratedIdentifierValueHandle(id, WomArrayType(womType))
        case GeneratedCallOutputValueHandle(first, second, womType) => GeneratedCallOutputValueHandle(first, second, WomArrayType(womType))
      } }
    }
  }

  implicit val IfElementUnlinkedValueGenerator: UnlinkedValueGenerator[IfElement] = new UnlinkedValueGenerator[IfElement] {
    override def generatedValueHandles(a: IfElement, typeAliases: Map[String, WomType], callables: Map[String, Callable]): ErrorOr[Set[GeneratedValueHandle]] = {
      a.graphElements.toList.traverse(_.generatedValueHandles(typeAliases, callables)).map(_.toSet.flatten) map { _.map {
        case GeneratedIdentifierValueHandle(id, womType) => GeneratedIdentifierValueHandle(id, WomOptionalType(womType).flatOptionalType)
        case GeneratedCallOutputValueHandle(first, second, womType) => GeneratedCallOutputValueHandle(first, second, WomOptionalType(womType).flatOptionalType)
      } }
    }
  }

  implicit val callElementUnlinkedValueGenerator: UnlinkedValueGenerator[CallElement] = new UnlinkedValueGenerator[CallElement] {
    override def generatedValueHandles(a: CallElement, typeAliases: Map[String, WomType], callables: Map[String, Callable]): ErrorOr[Set[GeneratedValueHandle]] = {
      def callableOutputToHandle(callable: Callable)(callableOutput: OutputDefinition): GeneratedValueHandle = {
        val callAlias = a.alias.getOrElse(callable.name)
        GeneratedCallOutputValueHandle(callAlias, callableOutput.name, callableOutput.womType)
      }

      callables.get(a.callableReference) match {
        case Some(callable) => callable.outputs.map(callableOutputToHandle(callable)).toSet.validNel
        case None => s"Cannot generate outputs for 'call ${a.callableReference}'. No such callable exists in [${callables.keySet.mkString(", ")}]".invalidNel
      }
    }
  }

  implicit val graphElementUnlinkedValueConsumer: GraphElementValueConsumer[WorkflowGraphElement] = new GraphElementValueConsumer[WorkflowGraphElement] {
    override def graphElementConsumedValueHooks(a: WorkflowGraphElement, typeAliases: Map[String, WomType], callables: Map[String, Callable]): ErrorOr[Set[UnlinkedConsumedValueHook]] = a match {
      case InputDeclarationElement(_, _, None) => Set.empty[UnlinkedConsumedValueHook].validNel
      case DeclarationElement(_, _, Some(expr)) => expr.expressionConsumedValueHooks.validNel
      case a: ScatterElement => a.graphElementConsumedValueHooks(typeAliases, callables)
      case a: IfElement => a.graphElementConsumedValueHooks(typeAliases, callables)
      case a: CallElement => a.graphElementConsumedValueHooks(typeAliases, callables)
      // TODO fill in other expression types
      case other => throw new Exception(s"Cannot generate consumed values for WorkflowGraphNodeElement $other")
    }
  }

  implicit val callElementUnlinkedValueConsumer: GraphElementValueConsumer[CallElement] = new GraphElementValueConsumer[CallElement] {
    override def graphElementConsumedValueHooks(a: CallElement, typeAliases: Map[String, WomType], callables: Map[String, Callable]): ErrorOr[Set[UnlinkedConsumedValueHook]] = {
      a.body match {
        case Some(callBodyElement: CallBodyElement) => callBodyElement.inputs.flatMap(_.expressionConsumedValueHooks).toSet.validNel
        case None => Set.empty[UnlinkedConsumedValueHook].valid
      }
    }
  }

  implicit val scatterElementUnlinkedValueConsumer: GraphElementValueConsumer[ScatterElement] = new GraphElementValueConsumer[ScatterElement] {
    override def graphElementConsumedValueHooks(a: ScatterElement, typeAliases: Map[String, WomType], callables: Map[String, Callable]): ErrorOr[Set[UnlinkedConsumedValueHook]] = {
      val bodyConsumedValuesValidation: ErrorOr[Set[UnlinkedConsumedValueHook]] = a.graphElements.toList.traverse(_.graphElementConsumedValueHooks(typeAliases, callables)).map(_.toSet.flatten)
      val scatterExpressionHooks: Set[UnlinkedConsumedValueHook] = a.scatterExpression.expressionConsumedValueHooks

      val bodyGeneratedValuesValidation: ErrorOr[Set[String]] = a.graphElements.toList.traverse(_.generatedValueHandles(typeAliases, callables)).map(_.toSet.flatten.map(_.linkableName))

      (bodyConsumedValuesValidation, bodyGeneratedValuesValidation) mapN { (bodyConsumedValues, bodyGeneratedValues) =>
        val unsatisfiedBodyElementHooks = bodyConsumedValues.filterNot {
          case UnlinkedIdentifierHook(id) => bodyGeneratedValues.contains(id) || id == a.scatterVariableName
          case UnlinkedCallOutputOrIdentifierAndMemberAccessHook(first, second) =>
            bodyGeneratedValues.contains(first) || bodyGeneratedValues.contains(s"$first.$second") || a.scatterVariableName == first
        }

        unsatisfiedBodyElementHooks ++ scatterExpressionHooks
      }
    }
  }

  implicit val ifElementUnlinkedValueConsumer: GraphElementValueConsumer[IfElement] = new GraphElementValueConsumer[IfElement] {
    override def graphElementConsumedValueHooks(a: IfElement, typeAliases: Map[String, WomType], callables: Map[String, Callable]): ErrorOr[Set[UnlinkedConsumedValueHook]] = {
      val bodyConsumedValuesValidation: ErrorOr[Set[UnlinkedConsumedValueHook]] = a.graphElements.toList.traverse(_.graphElementConsumedValueHooks(typeAliases, callables)).map(_.toSet.flatten)
      val scatterExpressionHooks: Set[UnlinkedConsumedValueHook] = a.conditionExpression.expressionConsumedValueHooks

      val bodyGeneratedValuesValidation: ErrorOr[Set[String]] = a.graphElements.toList.traverse(_.generatedValueHandles(typeAliases, callables)).map(_.toSet.flatten.map(_.linkableName))

      (bodyConsumedValuesValidation, bodyGeneratedValuesValidation) mapN { (bodyConsumedValues, bodyGeneratedValues) =>
        val unsatisfiedBodyElementHooks = bodyConsumedValues.filterNot {
          case UnlinkedIdentifierHook(id) => bodyGeneratedValues.contains(id)
          case UnlinkedCallOutputOrIdentifierAndMemberAccessHook(first, second) => bodyGeneratedValues.contains(first) || bodyGeneratedValues.contains(s"$first.$second")
        }

        unsatisfiedBodyElementHooks ++ scatterExpressionHooks
      }
    }
  }
}
