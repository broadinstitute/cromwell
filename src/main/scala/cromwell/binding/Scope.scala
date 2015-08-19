package cromwell.binding

import cromwell.binding.AstTools.AstNodeName
import cromwell.parser.WdlParser.{Ast, AstList}

import scala.annotation.tailrec
import scala.collection.mutable.{Map, HashMap}
import scala.collection.JavaConverters._

object Scope {

  /**
   * Generate a Scope Tree from an Ast.
   * Ast must match a Scope (currently: Workflow, Scatter or Call)
   * @throws UnsupportedOperationException if node does not match any supported Scopes
   * @return generated Scope Tree.
   */
  def generateTree(node: Ast, namespaces: Seq[WdlNamespace], tasks: Seq[Task], wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Scope = {

    val scopeIndexes: Map[Class[_ <: Scope], Int] = HashMap.empty

    /**
     * Retrieve the list of children Asts from the AST body.
     * Return empty seq if body is null or is not a list.
     */
    def getChildrenList(ast: Ast): Seq[Ast] = Option(ast.getAttribute("body")) collect {
      case list: AstList => list.asScala.toVector
    } getOrElse Seq.empty collect {
      case ast: Ast => ast
    }

    /**
     * Generate a scope from the ast parameter and all its descendants recursively.
     * @return Some(instance of Scope) if the node has a Scope equivalent, None otherwise
     */
    def generateScopeTree(node: Ast, namespaces: Seq[WdlNamespace], tasks: Seq[Task], wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter, parent: Option[Scope]): Option[Scope] = {
      val scope: Option[Scope] = node match {
        case w: Ast if w.getName == AstNodeName.Workflow => Option(Workflow(w, wdlSyntaxErrorFormatter, parent))
        case c: Ast if c.getName == AstNodeName.Call => Option(Call(c, namespaces, tasks, wdlSyntaxErrorFormatter, parent))
        case s: Ast if s.getName == AstNodeName.Scatter => Option(Scatter(s, scopeIndexes.getOrElse(classOf[Scatter], 0), parent))

        /* Cases need to be added here to support new scopes */
        case _ => None
      }
      //Generate and set children recursively
      scope foreach { sc =>
        //Increment index if needed
        scopeIndexes.put(sc.getClass, scopeIndexes.getOrElse(sc.getClass, 0) + 1)
        sc.setChildren(getChildrenList(node) flatMap { generateScopeTree(_, namespaces, tasks, wdlSyntaxErrorFormatter, Option(sc)) })
      }
      scope
    }

    generateScopeTree(node, namespaces, tasks, wdlSyntaxErrorFormatter, None)
      .getOrElse(throw new UnsupportedOperationException("Input node cannot be instantiated as a Scope"))
  }

  /**
   * Generate a workflow from an Ast. The ast must be one of a workflow.
   * @throws UnsupportedOperationException if the ast is not a workflow ast
   * @return a workflow
   */
  def generateWorkflow(ast: Ast, namespaces: Seq[WdlNamespace], tasks: Seq[Task], wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Workflow = {
    ast.getName match {
      case AstNodeName.Workflow => generateTree(ast, namespaces, tasks, wdlSyntaxErrorFormatter).asInstanceOf[Workflow]
      case nonWorkflowAst => throw new UnsupportedOperationException(s"Ast is not a 'Workflow Ast' but a '$nonWorkflowAst Ast'")
    }
  }

  /**
   * Collect Calls from a Seq of Scopes.
   * @param scopes scopes to loop through
   * @return Scopes instances that are Calls
   */
  def collectCalls(scopes: Seq[Scope]): Seq[Call] = scopes collect { case s: Call => s }

  /**
   * Collect all Calls from the given scope.
   * @param scopes scope to gather Calls from
   * @param calls for recursivity. Should be passed Nil in most cases.
   * @return all Calls inside the scope
   */
  @tailrec
  def collectAllCalls(scopes: Seq[Scope], calls: Seq[Call]): Seq[Call] = scopes match {
    case Nil => calls
    case l => collectAllCalls(l.flatMap(_.children), calls ++ collectCalls(l))
  }

  /**
   * Collect Scatters from a Seq of Scopes.
   * @param scopes scopes to loop through
   * @return Scopes instances that are Scatters
   */
  def collectScatters(scopes: Seq[Scope]): Seq[Scatter] = scopes collect { case s: Scatter => s }

  /**
   * Collect all Scatters from the given scope.
   * @param scopes scope to gather Scatters from
   * @param scatters for recursivity. Should be passed Nil in most cases.
   * @return all Scatters inside the scope
   */
  @tailrec
  def collectAllScatters(scopes: Seq[Scope], scatters: Seq[Scatter]): Seq[Scatter] = scopes match {
    case Nil => scatters
    case l => collectAllScatters(l.flatMap(_.children), scatters ++ collectScatters(l))
  }

  @tailrec
  def fullyQualifiedNameBuilder(scope: Option[Scope], fqn: String, fullDisplay: Boolean): String = scope match {
    case None => fqn.tail //Strip away the first "." of the name
    case Some(x: Scope) => fullyQualifiedNameBuilder(x.parent, (if (fullDisplay || x.appearsInFQN) s".${x.name}" else "") + fqn, fullDisplay)
  }

}

// FIXME: It'd be nice to have a notion of a parented and a not-parented Scope, see FIXME in Call about this
trait Scope {

  def name: String

  def appearsInFQN: Boolean = true

  val parent: Option[Scope]

  private var _children: Seq[Scope] = Seq.empty

  def children: Seq[Scope] = _children

  def fullyQualifiedName = Scope.fullyQualifiedNameBuilder(Option(this), "", fullDisplay = false)

  def fullyQualifiedNameWithIndexScopes = Scope.fullyQualifiedNameBuilder(Option(this), "", fullDisplay = true)

  def setChildren(children: Seq[Scope]) = {
    if (this._children.isEmpty) {
      this._children = children
    } else throw new UnsupportedOperationException("children is write-once")
  }

  /**
   * Convenience method to collect Calls from within a scope.
   * @return all calls contained in this scope (recursively)
   */
  def collectAllCalls = Scope.collectAllCalls(Seq(this), Nil)

  /**
   * Convenience method to collect Scatters from within a scope.
   * @return all scatters contained in this scope (recursively)
   */
  def collectAllScatters = Scope.collectAllScatters(Seq(this), Nil)

}
