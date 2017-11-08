package wdl

import cats.data.Validated.Valid
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl.AstTools.EnhancedAstNode
import wdl4s.parser.WdlParser.{Ast, AstNode}
import wom.WorkflowInput
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
    case None => Option(TaskInput(unqualifiedName, womType))
  }

  def asWorkflowInput: Option[WorkflowInput] = expression match {
    case Some(_) => None
    case None => Some(WorkflowInput(fullyQualifiedName, womType))
  }

  def toWdlString: String = {
    val expr = expression.map(e => s" = ${e.toWomString}").getOrElse("")
    s"${womType.toDisplayString} $unqualifiedName$expr"
  }

  final lazy val upstreamReferences = expression.toSeq.flatMap(_.variableReferences)

  override def toString: String = {
    s"[Declaration type=${womType.toDisplayString} name=$unqualifiedName expr=${expression.map(_.toWomString)}]"
  }
}

object Declaration {

  sealed trait WdlDeclarationNode {
    lazy val toGraphNode: GraphNode = this match {
      case InputDeclarationNode(graphInputNode) => graphInputNode
      case IntermediateValueDeclarationNode(expressionNode) => expressionNode
      case GraphOutputDeclarationNode(graphOutputNode) => graphOutputNode
    }
    lazy val singleOutputPort: Option[GraphNodePort.OutputPort] = this match {
      case InputDeclarationNode(graphInputNode) => Option(graphInputNode.singleOutputPort)
      case IntermediateValueDeclarationNode(expressionNode) => Option(expressionNode.singleExpressionOutputPort)
      case GraphOutputDeclarationNode(_) => None
    }
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

  def buildWomNode(decl: DeclarationInterface, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[WdlDeclarationNode] = {

    def declarationAsExpressionNode(wdlExpression: WdlExpression) = {
      val womExpression = WdlWomExpression(wdlExpression, None)
      for {
        inputMapping <- WdlWomExpression.findInputsforExpression(womExpression, localLookup, outerLookup, preserveIndexForOuterLookups)
        expressionNode <- ExposedExpressionNode.fromInputMapping(decl.womIdentifier, womExpression, decl.womType, inputMapping)
      } yield IntermediateValueDeclarationNode(expressionNode)
    }

    def workflowOutputAsGraphOutputNode(wdlExpression: WdlExpression) = {
      val womExpression = WdlWomExpression(wdlExpression, None)
      for {
        inputMapping <- WdlWomExpression.findInputsforExpression(womExpression, localLookup, outerLookup, preserveIndexForOuterLookups)
        graphOutputNode <- ExpressionBasedGraphOutputNode.fromInputMapping(decl.womIdentifier, womExpression, decl.womType, inputMapping)
      } yield GraphOutputDeclarationNode(graphOutputNode)
    }

    decl match {
      case Declaration(opt: WomOptionalType, _, None, _, _) => Valid(InputDeclarationNode(OptionalGraphInputNode(decl.womIdentifier, opt)))
      case Declaration(_, _, None, _, _) => Valid(InputDeclarationNode(RequiredGraphInputNode(decl.womIdentifier, decl.womType)))
      case Declaration(_, _, Some(expr), _, _) => declarationAsExpressionNode(expr)
      case WorkflowOutput(_, _, expr, _, _) => workflowOutputAsGraphOutputNode(expr)
    }
  }
}

case class Declaration(womType: WomType,
                       unqualifiedName: String,
                       expression: Option[WdlExpression],
                       override val parent: Option[Scope],
                       ast: Ast) extends DeclarationInterface
