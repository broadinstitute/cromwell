package wdl.draft2.model

import cats.data.Validated.Valid
import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl.draft2.model
import wdl.draft2.model.AstTools.EnhancedAstNode
import wdl.draft2.parser.WdlParser.{Ast, AstNode}
import wom.callable.Callable.{InputDefinition, OverridableInputDefinitionWithDefault, OptionalInputDefinition, RequiredInputDefinition}
import wom.graph._
import wom.graph.expression.{ExposedExpressionNode, ExpressionNode}
import wom.types.{WomArrayType, WomOptionalType, WomType}

object DeclarationInterface {
  /**
    * Depending on who is asking, the type of a declaration can vary.
    * e.g
    * task a {
    *   ...
    *   output {
    *     String o
    *   }
    * }
    * workflow w {
    *   scatter (...) {
    *     call a
    *     String s = a.o # Inside the scatter a.o is a String
    *   }
    *
    *   Array[String] s2 = a.o # Outside the scatter it's an Array[String]
    * }
    */
  def relativeWdlType(from: Scope, target: DeclarationInterface, womType: WomType): WomType = {
    target.closestCommonAncestor(from) map { ancestor =>
      target.ancestrySafe.takeWhile(_ != ancestor).foldLeft(womType){
        case (acc, _: Scatter) => WomArrayType(acc)
        case (acc, _: If) => WomOptionalType(acc).flatOptionalType
        case (acc, _) => acc
      }
    } getOrElse womType
  }
}

/**
  * Represents a declaration which can show up in a workflow or a task context.  For example
  *
  * task test {
  *   File test_file
  *   command { ... }
  * }
  *
  * workflow wf {
  *   String wf_string = "abc"
  *   call test { input: s=wf_string }
  * }
  *
  * Both the definition of test_file and wf_string are declarations
  */
trait DeclarationInterface extends WdlGraphNodeWithUpstreamReferences {
  def womType: WomType
  def expression: Option[WdlExpression]
  def ast: Ast

  def relativeWdlType(from: Scope): WomType = DeclarationInterface.relativeWdlType(from, this, womType)

  def asTaskInput: Option[TaskInput] = expression match {
    case Some(_) => None
    case None => Option(model.TaskInput(unqualifiedName, womType))
  }

  /**
    * Decide whether this declaration should be exposed as a workflow input, and if so, what kind
    */
  def asWorkflowInput: Option[InputDefinition] = (expression, womType) match {
    case (None, optionalType: WomOptionalType) =>
      Some(OptionalInputDefinition(fullyQualifiedName, optionalType.flatOptionalType))
    case (None, nonOptionalType) =>
      Some(RequiredInputDefinition(fullyQualifiedName, nonOptionalType))
    // We only make declarations with expressions into inputs if they don't depend on previous tasks or decls:
    case (Some(expr), other) if upstreamAncestry.isEmpty =>
      Some(OverridableInputDefinitionWithDefault(fullyQualifiedName, other, WdlWomExpression(expr, this)))
    case _ =>
      None
  }

  def toWdlString: String = {
    val expr = expression.map(e => s" = ${e.toWomString}").getOrElse("")
    s"${womType.stableName} $unqualifiedName$expr"
  }

  final lazy val upstreamReferences = expression.toSeq.flatMap(_.variableReferences(this))

  override def toString: String = {
    s"[Declaration type=${womType.stableName} name=$unqualifiedName expr=${expression.map(_.toWomString)}]"
  }
}

object Declaration {

  sealed trait WdlDeclarationNode {
    lazy val toGraphNode: GraphNode = this match {
      case InputDeclarationNode(graphInputNode) => graphInputNode
      case IntermediateValueDeclarationNode(expressionNode) => expressionNode
      case GraphOutputDeclarationNode(graphOutputNode) => graphOutputNode
    }
    lazy val singleOutputPort = this match {
      case InputDeclarationNode(graphInputNode) => graphInputNode.singleOutputPort
      case IntermediateValueDeclarationNode(expressionNode) => expressionNode.singleOutputPort
      case GraphOutputDeclarationNode(outputNode) => outputNode.graphOutputPort
    }
    lazy val localName = toGraphNode.localName
  }
  final case class InputDeclarationNode(graphInputNode: GraphInputNode) extends WdlDeclarationNode
  final case class IntermediateValueDeclarationNode(expressionNode: ExpressionNode) extends WdlDeclarationNode
  final case class GraphOutputDeclarationNode(graphOutputNode: GraphOutputNode) extends WdlDeclarationNode

  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): Declaration = {
    Declaration(
      ast.getAttribute("type").womType(wdlSyntaxErrorFormatter),
      ast.getAttribute("name").sourceString,
      ast.getAttribute("expression") match {
        case a: AstNode => Option(WdlExpression(a))
        case _ => None
      },
      parent,
      ast
    )
  }

  def buildWdlDeclarationNode(decl: DeclarationInterface, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[WdlDeclarationNode] = {

    def declarationAsExpressionNode(wdlExpression: WdlExpression) = {
      val womExpression = WdlWomExpression(wdlExpression, decl)
      for {
        inputMapping <- WdlWomExpression.findInputsforExpression(womExpression, localLookup, outerLookup, preserveIndexForOuterLookups, decl)
        expressionNode <- ExposedExpressionNode.fromInputMapping(decl.womIdentifier, womExpression, decl.womType, inputMapping)
      } yield IntermediateValueDeclarationNode(expressionNode)
    }

    def workflowOutputAsGraphOutputNode(wdlExpression: WdlExpression) = {
      val womExpression = WdlWomExpression(wdlExpression, decl)
      for {
        inputMapping <- WdlWomExpression.findInputsforExpression(womExpression, localLookup, outerLookup, preserveIndexForOuterLookups, decl)
        graphOutputNode <- ExpressionBasedGraphOutputNode.fromInputMapping(decl.womIdentifier, womExpression, decl.womType, inputMapping)
      } yield GraphOutputDeclarationNode(graphOutputNode)
    }

    def asWorkflowInput(inputDefinition: InputDefinition): GraphInputNode = inputDefinition match {
      case RequiredInputDefinition(_, womType, _, _) => RequiredGraphInputNode(decl.womIdentifier, womType, decl.womIdentifier.fullyQualifiedName.value)
      case OptionalInputDefinition(_, optionalType, _, _) => OptionalGraphInputNode(decl.womIdentifier, optionalType, decl.womIdentifier.fullyQualifiedName.value)
      case OverridableInputDefinitionWithDefault(_, womType, default, _, _) => OptionalGraphInputNodeWithDefault(decl.womIdentifier, womType, default, decl.womIdentifier.fullyQualifiedName.value)
      case other => throw new RuntimeException(s"Programmer Error! If you got here you probably changed draft 2 to try to do some draft 3 like things, but this draft 2 function isn't set up to produce or handle ${other.getClass.getSimpleName} yet!")
    }

    (decl.asWorkflowInput, decl) match {
      case (_, WorkflowOutput(_, _, expr, _, _)) => workflowOutputAsGraphOutputNode(expr)
      case (Some(inputDefinition), _) => Valid(InputDeclarationNode(asWorkflowInput(inputDefinition)))
      case (None, Declaration(_, _, Some(expr), _, _)) => declarationAsExpressionNode(expr)
      case (_, other) => s"Non-input declaration ($other) didn't have an expression".invalidNel
    }
  }
}

case class Declaration(womType: WomType,
                       unqualifiedName: String,
                       expression: Option[WdlExpression],
                       override val parent: Option[Scope],
                       ast: Ast) extends DeclarationInterface
