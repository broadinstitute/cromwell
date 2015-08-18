package cromwell.binding

import scala.annotation.tailrec
import scala.language.postfixOps

object Scope {
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
  def fullyQualifiedNameBuilder(scope: Option[Scope], fqn: String, fullDisplay: Boolean, leaf: Boolean): String = {
    scope match {
      case Some(x: Scope) =>
        fullyQualifiedNameBuilder(
          x.parent,
          (if (fullDisplay || x.appearsInFqn || leaf) s".${x.name}" else "") + fqn,
          fullDisplay,
          leaf = false)
      case None => fqn.tail //Strip away the first "." of the name
    }
  }

}

trait Scope {
  def name: LocallyQualifiedName
  def appearsInFqn: Boolean = true

  val parent: Option[Scope]
  private var _children: Seq[Scope] = Seq.empty
  def children: Seq[Scope] = _children

  def children_=[Child <: Scope](children: Seq[Child]): Unit = {
    if (this._children.isEmpty) {
      this._children = children
    } else throw new UnsupportedOperationException("children is write-once")
  }

  def fullyQualifiedName =
    Scope.fullyQualifiedNameBuilder(Option(this), "", fullDisplay = false, leaf = true)

  def fullyQualifiedNameWithIndexScopes =
    Scope.fullyQualifiedNameBuilder(Option(this), "", fullDisplay = true, leaf = true)

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

  /*
   * Calls and scatters are accessed frequently so this avoids traversing the whole children tree every time.
   * Lazy because children are not provided at instantiation but rather later during tree building process.
   * This prevents evaluation from being done before children have been set.
   *
   * FIXME: In a world where Scope wasn't monolithic, these would be moved around
   */
  lazy val calls: Seq[Call] = collectAllCalls
  lazy val scatters: Seq[Scatter] = collectAllScatters

  /**
   * Recurses up the scope tree until it hits one w/o a parent.
   *
   * FIXME: In a world where Scope wasn't monolithic this would traverse until it hit a RootScope or whatever
   */
  final def rootScope: Scope = this.parent map { _.rootScope } getOrElse this

  // FIXME: In a world where Scope wasn't monolithic, these would be moved out of here
  def prerequisiteScopes: Set[Scope]
  def prerequisiteCallNames: Set[LocallyQualifiedName]

  /**
   *  Returns a set of Calls corresponding to the prerequisiteCallNames
   *
   *  Dropping any unfound Calls to the floor but we're already validating that all calls are sane at ingest.
   *  It's icky because it relies on that validation not changing, but ...
   */
  lazy val prerequisiteCalls: Set[Scope] = prerequisiteCallNames flatMap rootScope.callByName
  def callByName(callName: LocallyQualifiedName): Option[Call] = calls find { _.name == callName }

  def ancestry: Seq[Scope] = parent match {
    case Some(p) => Seq(p) ++ p.ancestry
    case None => Seq.empty[Scope]
  }

  def closestCommonAncestor(other: Scope): Option[Scope] = {
    val otherAncestry = other.ancestry
    ancestry find { otherAncestry.contains(_) }
  }
}
