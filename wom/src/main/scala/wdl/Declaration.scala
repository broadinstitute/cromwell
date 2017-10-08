package wdl

import cats.data.Validated.Valid
import lenthall.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import wdl4s.parser.WdlParser.{Ast, AstNode}
import wdl.AstTools.EnhancedAstNode
import wdl.types.{WdlArrayType, WdlOptionalType, WdlType}
import wom.graph._

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
  def relativeWdlType(from: Scope, target: DeclarationInterface, wdlType: WdlType): WdlType = {
    target.closestCommonAncestor(from) map { ancestor =>
      target.ancestrySafe.takeWhile(_ != ancestor).foldLeft(wdlType){
        case (acc, _: Scatter) => WdlArrayType(acc)
        case (acc, _: If) => WdlOptionalType(acc)
        case (acc, _) => acc
      }
    } getOrElse wdlType
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
  def wdlType: WdlType
  def expression: Option[WdlExpression]
  def ast: Ast

  def relativeWdlType(from: Scope): WdlType = DeclarationInterface.relativeWdlType(from, this, wdlType)

  def asTaskInput: Option[TaskInput] = expression match {
    case Some(_) => None
    case None => Option(TaskInput(unqualifiedName, wdlType))
  }

  def asWorkflowInput: Option[WorkflowInput] = expression match {
    case Some(_) => None
    case None => Some(WorkflowInput(fullyQualifiedName, wdlType))
  }

  def toWdlString: String = {
    val expr = expression.map(e => s" = ${e.toWdlString}").getOrElse("")
    s"${wdlType.toWdlString} $unqualifiedName$expr"
  }

  final lazy val upstreamReferences = expression.toSeq.flatMap(_.variableReferences)

  override def toString: String = {
    s"[Declaration type=${wdlType.toWdlString} name=$unqualifiedName expr=${expression.map(_.toWdlString)}]"
  }
}

object Declaration {

  sealed trait DeclarationNode {
    def toGraphNode = this match {
      case InputDeclarationNode(graphInputNode) => graphInputNode
      case IntermediateValueDeclarationNode(expressionNode) => expressionNode
      case GraphOutputDeclarationNode(graphOutputNode) => graphOutputNode
    }
    def singleOutputPort: Option[GraphNodePort.OutputPort] = this match {
      case InputDeclarationNode(graphInputNode) => Option(graphInputNode.singleOutputPort)
      case IntermediateValueDeclarationNode(expressionNode) => Option(expressionNode.singleExpressionOutputPort)
      case GraphOutputDeclarationNode(_) => None
    }
  }
  final case class InputDeclarationNode(graphInputNode: GraphInputNode) extends DeclarationNode
  final case class IntermediateValueDeclarationNode(expressionNode: ExpressionNode) extends DeclarationNode
  final case class GraphOutputDeclarationNode(graphOutputNode: GraphOutputNode) extends DeclarationNode

  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): Declaration = {
    Declaration(
      ast.getAttribute("type").wdlType(wdlSyntaxErrorFormatter),
      ast.getAttribute("name").sourceString,
      ast.getAttribute("expression") match {
        case a: AstNode => Option(WdlExpression(a))
        case _ => None
      },
      parent,
      ast
    )
  }

  def buildWomNode(decl: DeclarationInterface, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, GraphNodePort.OutputPort]): ErrorOr[DeclarationNode] = {

    val inputName = decl.unqualifiedName

    def declarationAsExpressionNode(wdlExpression: WdlExpression) = {
      val womExpression = WdlWomExpression(wdlExpression, None)
      for {
        uninstantiatedExpression <- WdlWomExpression.findInputsforExpression(inputName, womExpression, localLookup, outerLookup)
        expressionNode <- ExpressionNode.linkWithInputs(decl.womIdentifier, womExpression, uninstantiatedExpression.inputMapping)
      } yield IntermediateValueDeclarationNode(expressionNode)
    }

    def workflowOutputAsGraphOutputNode(wdlExpression: WdlExpression) = {
      val womExpression = WdlWomExpression(wdlExpression, None)
      for {
        uninstantiatedExpression <- WdlWomExpression.findInputsforExpression(inputName, womExpression, localLookup, outerLookup)
        graphOutputNode <- ExpressionBasedGraphOutputNode.linkWithInputs(decl.womIdentifier, decl.wdlType, womExpression, uninstantiatedExpression.inputMapping)
      } yield GraphOutputDeclarationNode(graphOutputNode)
    }

    decl match {
      case Declaration(opt: WdlOptionalType, _, None, _, _) => Valid(InputDeclarationNode(OptionalGraphInputNode(decl.womIdentifier, opt)))
      case Declaration(_, _, None, _, _) => Valid(InputDeclarationNode(RequiredGraphInputNode(decl.womIdentifier, decl.wdlType)))
      case Declaration(_, _, Some(expr), _, _) => declarationAsExpressionNode(expr)
      case WorkflowOutput(_, _, expr, _, _) => workflowOutputAsGraphOutputNode(expr)
    }
  }
}

case class Declaration(wdlType: WdlType,
                       unqualifiedName: String,
                       expression: Option[WdlExpression],
                       override val parent: Option[Scope],
                       ast: Ast) extends DeclarationInterface
