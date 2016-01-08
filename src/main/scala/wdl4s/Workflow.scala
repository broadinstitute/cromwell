package wdl4s

import wdl4s.AstTools.{AstNodeName, EnhancedAstNode}
import wdl4s.types.WdlType
import wdl4s.parser.WdlParser.{Ast, AstList, SyntaxError, Terminal}

import scala.collection.JavaConverters._
import scala.collection.mutable
import scala.language.postfixOps

object Workflow {

  def apply(ast: Ast,
            namespaces: Seq[WdlNamespace],
            tasks: Seq[Task],
            wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Workflow = {
    ast.getName match {
      case AstNodeName.Workflow =>
        generateWorkflowScopes(ast, namespaces, tasks, wdlSyntaxErrorFormatter)
      case nonWorkflowAst =>
        throw new UnsupportedOperationException(s"Ast is not a 'Workflow Ast' but a '$nonWorkflowAst Ast'")
    }
  }

  private def generateWorkflowScopes(workflowAst: Ast,
                                     namespaces: Seq[WdlNamespace],
                                     tasks: Seq[Task],
                                     wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Workflow = {

    /**
     * Incremented indexes per type of class. Currently only used by Scatter.
     * Indexes are incremented just before use + recursion (prefix, vs. infix or postfix),
     * thus "start" at zero, but default to -1.
     */
    val scopeIndexes: mutable.Map[Class[_ <: Scope], Int] = mutable.HashMap.empty.withDefaultValue(-1)

    /**
     * Retrieve the list of children ASTs from the AST body.
     * Return empty seq if body is null or is not a list.
     */
    def getChildrenList(ast: Ast): Seq[Ast] = {
      val body = Option(ast.getAttribute("body")).filter(_.isInstanceOf[AstList]).map(_.astListAsVector)
      body match {
        case Some(asts) => asts collect { case node: Ast if Seq(AstNodeName.Call, AstNodeName.Scatter).contains(node.getName) => node }
        case None => Seq.empty[Ast]
      }
    }

    /**
     * Sets the child scopes on this parent from the ast.
     * @param parentAst Parent ast.
     * @param parentScope Parent scope.
     */
    def setChildScopes(parentAst: Ast, parentScope: Scope): Unit = {
      parentScope.children = for {
        child <- getChildrenList(parentAst)
        scope = generateScopeTree(child, namespaces, tasks, wdlSyntaxErrorFormatter, parentScope)
      } yield scope
    }

    /**
     * Generate a scope from the ast parameter and all its descendants recursively.
     * @return Scope equivalent
     */
    def generateScopeTree(node: Ast,
                          namespaces: Seq[WdlNamespace],
                          tasks: Seq[Task],
                          wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter,
                          parent: Scope): Scope = {

      val scope: Scope = node.getName match {
        case AstNodeName.Call =>
          Call(node, namespaces, tasks, wdlSyntaxErrorFormatter, Option(parent))
        case AstNodeName.Scatter =>
          scopeIndexes(classOf[Scatter]) += 1
          Scatter(node, scopeIndexes(classOf[Scatter]), Option(parent))
      }
      // Generate and set children recursively
      setChildScopes(node, scope)
      scope
    }

    val workflow = Workflow(workflowAst, wdlSyntaxErrorFormatter)
    setChildScopes(workflowAst, workflow)
    workflow
  }

  private def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Workflow = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val declarations = ast.findAsts(AstNodeName.Declaration).map(Declaration(_, name, wdlSyntaxErrorFormatter))
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

    new Workflow(name, declarations, workflowOutputsDecls)
  }
}

/**
 * Represents a `workflow` definition in WDL which currently
 * only supports a set of `call` declarations and a name for
 * the workflow
 *
 * @param unqualifiedName The name of the workflow
 */
case class Workflow(unqualifiedName: String,
                    declarations: Seq[Declaration],
                    workflowOutputDecls: Seq[WorkflowOutputDeclaration]) extends Executable with Scope {
  override val prerequisiteScopes = Set.empty[Scope]
  override val prerequisiteCallNames = Set.empty[String]
  override val parent: Option[Scope] = None

  /**
   * FQNs for all inputs to this workflow and their associated types and possible postfix quantifiers.
   *
   * @return a Map[FullyQualifiedName, WorkflowInput] representing the
   *         inputs that the user needs to provide to this workflow
   */
  def inputs: Map[FullyQualifiedName, WorkflowInput] = {
    val callInputs = for { call <- calls; input <- call.unsatisfiedInputs } yield input
    val declarationInputs = for { declaration <- declarations; input <- declaration.asWorkflowInput } yield input
    (callInputs ++ declarationInputs) map { input => input.fqn -> input } toMap
  }

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
      call: Call <- calls
      output <- call.task.outputs
    } yield PotentialReportableSymbol(s"${call.fullyQualifiedName}.${output.name}", output.wdlType, matchWorkflowOutputWildcards = true)

    val inputs: Seq[PotentialReportableSymbol] = for {
      call: Call <- calls
      input <- call.task.inputs
    } yield PotentialReportableSymbol(s"${call.fullyQualifiedName}.${input.name}", input.wdlType, matchWorkflowOutputWildcards = false)

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

  override def rootWorkflow: Workflow = this
}

case class ReportableSymbol(fullyQualifiedName: FullyQualifiedName, wdlType: WdlType)