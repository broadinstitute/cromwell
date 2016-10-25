package wdl4s

import wdl4s.AstTools.{AstNodeName, EnhancedAstNode}
import wdl4s.types.WdlType
import wdl4s.parser.WdlParser.{Ast, AstList, SyntaxError, Terminal}

import scala.collection.JavaConverters._
import scala.collection.mutable
import scala.language.postfixOps

object Workflow {
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Workflow = {
    if (ast.getName != AstNodeName.Workflow) {
      throw new UnsupportedOperationException(s"Expecting Workflow AST, got a ${ast.getName} AST")
    }
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val callNames = ast.findAsts(AstNodeName.Call).map {call =>
      Option(call.getAttribute("alias")).getOrElse(call.getAttribute("task"))
    }
    val workflowOutputsDecls = ast.findAsts(AstNodeName.WorkflowOutput) map { wfOutput =>
      val wildcard = Option(wfOutput.getAttribute("wildcard")).map(_.sourceString).getOrElse("").nonEmpty
      val outputFqn = name + "." + wfOutput.getAttribute("fqn").sourceString
      WorkflowOutputDeclaration(outputFqn, wildcard)
    }

    callNames groupBy { _.sourceString } foreach {
      case (_, terminals) if terminals.size > 1 =>
        throw new SyntaxError(wdlSyntaxErrorFormatter.multipleCallsAndHaveSameName(terminals.asInstanceOf[Seq[Terminal]]))
      case _ =>
    }

    new Workflow(name, workflowOutputsDecls, ast)
  }
}

case class Workflow(unqualifiedName: String,
                    workflowOutputDecls: Seq[WorkflowOutputDeclaration],
                    ast: Ast) extends Scope {

  /**
   * FQNs for all inputs to this workflow and their associated types and possible postfix quantifiers.
   *
   * @return a Map[FullyQualifiedName, WorkflowInput] representing the
   *         inputs that the user needs to provide to this workflow
   */
  def inputs: Map[FullyQualifiedName, WorkflowInput] = {
    val callInputs = for {
      call <- calls
      input <- call.unsatisfiedInputs
    } yield input

    val declarationInputs = for {
      declaration <- declarations
      input <- declaration.asWorkflowInput
    } yield input

    (callInputs ++ declarationInputs) map { input => input.fqn -> input } toMap
  }

  /** First tries to find any Call with name `name`.  If not found,
    * Fallback to looking at immediate children or delegating to parent node
    */
  override def resolveVariable(name: String, relativeTo: Scope = this): Option[Scope with GraphNode] = {
    findCallByName(name) match {
      case call: Some[Call] => call
      case _ => super.resolveVariable(name, relativeTo)
    }
  }

  /**
    * Find a Call object by name.  For example,
    *
    * {{{
    * workflow w {
    *   call foobar
    * }
    * }}}
    *
    * calling findCallByName("foobar") will return a Some(call).  All
    * other values would return a None
    *
    * @param name name of call to return
    * @return Some(Call) if one with that name was found otherwise None
    */
  def findCallByName(name: String): Option[Call] = calls.find(_.unqualifiedName == name)

  /**
   * All outputs for this workflow and their associated types
   *
   * @return a Map[FullyQualifiedName, WdlType] representing the union
   *         of all outputs from all `call`s within this workflow
   */
  lazy val outputs: Seq[ReportableSymbol] = {

    case class PotentialReportableSymbol(name: String, wdlType: WdlType, matchWorkflowOutputWildcards: Boolean)

    // Build a list of ALL potentially reportable symbols, and whether they're allowed to match
    // wildcards in the workflow's output {...} spec.
    val outputs: Seq[PotentialReportableSymbol] = for {
      call: Call <- calls.toSeq
      output <- call.task.outputs
    } yield PotentialReportableSymbol(s"${call.fullyQualifiedName}.${output.unqualifiedName}", output.wdlType, matchWorkflowOutputWildcards = true)

    val inputs: Seq[PotentialReportableSymbol] = for {
      call: Call <- calls.toSeq
      input <- call.task.declarations
    } yield PotentialReportableSymbol(s"${call.fullyQualifiedName}.${input.unqualifiedName}", input.wdlType, matchWorkflowOutputWildcards = false)

    val filtered = if (workflowOutputDecls isEmpty) {
      outputs
    } else {
      (outputs ++ inputs) filter {
        case PotentialReportableSymbol(fqn, wdlType, wildcardsAllowed) =>
          workflowOutputDecls.isEmpty || workflowOutputDecls.exists(_.outputMatchesDeclaration(fqn, wildcardsAllowed))
      }
    }

    filtered map { case PotentialReportableSymbol(fqn, value, wildcardAllowed) => ReportableSymbol(fqn, value) }
  }

  override def toString = s"[Workflow $fullyQualifiedName]"
}

case class ReportableSymbol(fullyQualifiedName: FullyQualifiedName, wdlType: WdlType)