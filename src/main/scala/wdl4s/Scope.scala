package wdl4s

import wdl4s.expression.WdlFunctions
import wdl4s.parser.WdlParser.Ast
import wdl4s.values.{WdlArray, WdlValue}

import scala.language.postfixOps
import scala.util.{Failure, Success}

trait Scope {
  def unqualifiedName: LocallyQualifiedName
  def appearsInFqn: Boolean = true
  def ast: Ast

  /**
    * Parent scope
    */
  def parent: Option[Scope] = _parent
  private var _parent: Option[Scope] = None
  def parent_=[Child <: Scope](scope: Scope): Unit = {
    if (this._parent.isEmpty) this._parent = Option(scope)
    else throw new UnsupportedOperationException("parent is write-once")
  }

  /**
    * Child scopes, in the order that they appear in the source code
    */
  def children: Seq[Scope] = _children
  private var _children: Seq[Scope] = Seq.empty
  def children_=[Child <: Scope](children: Seq[Child]): Unit = {
    if (this._children.isEmpty) {
      this._children = children
    } else throw new UnsupportedOperationException("children is write-once")
  }

  /**
    * Containing namespace
    */
  def namespace: WdlNamespace = _namespace
  private var _namespace: WdlNamespace = null
  def namespace_=[Child <: WdlNamespace](ns: WdlNamespace): Unit = {
    if (Option(this._namespace).isEmpty) {
      this._namespace = ns
    } else throw new UnsupportedOperationException("namespace is write-once")
  }

  /**
    * Seq(parent, grandparent, great grandparent, ..., WdlNamespace)
    */
  lazy val ancestry: Seq[Scope] = parent match {
    case Some(p) => Seq(p) ++ p.ancestry
    case None => Seq.empty[Scope]
  }

  // This is needed in one specific case during WdlNamespace validation
  // where we need to compute the ancestries at a point where the full
  // parent branch has not been set yet.
  private [wdl4s] def ancestrySafe: Seq[Scope] = parent match {
    case Some(p) => Seq(p) ++ p.ancestrySafe
    case None => Seq.empty[Scope]
  }
  
  /**
    * All children ++ children's children ++ etc
    */
  lazy val descendants: Set[Scope] = (children ++ children.flatMap(_.descendants)).toSet

  /**
    * Descendants that are Calls
    */
  lazy val calls: Set[Call] = descendants.collect({ case c: Call => c })
  
  lazy val taskCalls: Set[TaskCall] = calls collect { case c: TaskCall => c }
  
  lazy val workflowCalls: Set[WorkflowCall] = calls collect { case c: WorkflowCall => c }

  /**
    * Descendants that are Scatters
    */
  lazy val scatters: Set[Scatter] = descendants.collect({ case s: Scatter => s })

  /**
    * Declarations within this Scope, in the order that they appear in source code
    */
  lazy val declarations: Seq[Declaration] = children.collect({ case d: Declaration => d})

  /**
    * String identifier for this scope.  this.namespace.resolve(this.fullyQualifiedName) == this
    */
  def fullyQualifiedName = {
    (ancestry.reverse.filter(_.appearsInFqn).map(_.unqualifiedName) :+ unqualifiedName).mkString(".")
  }

  /**
    * Similar to fullyQualifiedName but relatively to an ancestry scope.
    * e.g.
    * Workflow w
    *   Call a
    *     Output o
    *     
    * o.locallyQualified(a) = "a.o"
    * o.locallyQualified(w) = o.fullyQualifiedName = "w.a.o"
    */
  def locallyQualifiedName(relativeTo: Scope): String = {
    // Take ancestries until we reach relativeTo
    (ancestry.takeWhile(_ != relativeTo)
      // we want relativeTo in the lqn but it's been rejected by takeWhile wo add it back
      .:+(relativeTo)
      // Reverse because we start from the scope and climb up the ancestry tree but in the end we want a top-bottom lqn 
      .reverse
      // Get rid of scatters, ifs... because we don't want them here
      .filter(_.appearsInFqn)
      // Take the unqualifiedName of each scope
      .map(_.unqualifiedName)
      // Add the current scope
      :+ unqualifiedName)
      // Concatenate all of this
      .mkString(".")
  }
  
  /**
    * String identifier for this scope, with hidden scope information.
    *
    * this.namespace.resolve(this.fullyQualifiedNameWithIndexScopes) == this
    */
  def fullyQualifiedNameWithIndexScopes = {
    (Seq(this) ++ ancestry).reverse.map(_.unqualifiedName).filter(_.nonEmpty).mkString(".")
  }

  /**
    * Given another scope, returns the closest common ancestor between the two scopes,
    * if one exists at all
    *
    * @return closest common ancestor
    */
  def closestCommonAncestor(other: Scope): Option[Scope] = {
    val otherAncestry = other.ancestry
    ancestry find { otherAncestry.contains(_) }
  }

  /**
    * Performs scope resolution starting from this scope and walking up the lexical hierarchy
    * until it finds a GraphNode with the `name` as its unqualifiedName
    */
  def resolveVariable(name: String, relativeTo: Scope = this): Option[Scope with GraphNode] = {
    val siblingScopes = if (children.contains(relativeTo))
      // For declarations, only resolve to declarations that are lexically before this declaration
      children.dropRight(children.size - children.indexOf(relativeTo) )
    else children

    val localLookup = siblingScopes collect {
      case d: Declaration if d.unqualifiedName == name => d
      case c: TaskCall if c.unqualifiedName == name => c
      case o: TaskOutput if o.unqualifiedName == name => o
    }

    // If this is a scatter and the variable being resolved is the item
    val scatterLookup = Seq(this) collect {
      case s: Scatter if s.item == name => s
    }

    (scatterLookup ++ localLookup).headOption match {
      case scope: Some[_] => scope
      case None => parent.flatMap(_.resolveVariable(name, relativeTo))
    }
  }

  /**
    * This will return a lookup function for evaluating expressions which will traverse up the
    * scope hierarchy to find a value for `name`.  An exception will be thrown if a value cannot
    * be found for `name`
    *
    * @param knownInputs All known values of FQNs
    * @param wdlFunctions Implementation of WDL functions for expression evaluation
    * @param shards For resolving specific shards of scatter blocks
    * @return String => WdlValue lookup function rooted at `scope`
    * @throws VariableNotFoundException => If no errors occurred, but also `name` didn't resolve to any value
    * @throws VariableLookupException if anything else goes wrong in looking up a value for `name`
    */
  def lookupFunction(knownInputs: WorkflowCoercedInputs,
                     wdlFunctions: WdlFunctions[WdlValue],
                     outputResolver: OutputResolver = NoOutputResolver,
                     shards: Map[Scatter, Int] = Map.empty[Scatter, Int],
                     relativeTo: Scope = this): String => WdlValue = {

    def handleScatterResolution(scatter: Scatter): Option[WdlValue] = {
      // This case will happen if `name` references a Scatter.item (i.e. `x` in expression scatter(x in y) {...})
      val evaluatedCollection = scatter.collection.evaluate(scatter.lookupFunction(knownInputs, wdlFunctions, outputResolver, shards), wdlFunctions)
      val scatterShard = shards.get(scatter)

      (evaluatedCollection, scatterShard) match {
        case (Success(value: WdlArray), Some(shard)) if 0 <= shard && shard < value.value.size =>
          value.value.lift(shard)
        case (Success(value: WdlArray), Some(shard)) =>
          throw new VariableLookupException(s"Scatter expression (${scatter.collection.toWdlString}) evaluated to an array of ${value.value.size} elements, but element ${shard} was requested.")
        case (Success(value: WdlValue), _) =>
          throw new VariableLookupException(s"Expected scatter expression (${scatter.collection.toWdlString}) to evaluate to an Array.  Instead, got a ${value}")
        case (Failure(ex), _) =>
          throw new VariableLookupException(s"Failed to evaluate scatter expression (${scatter.collection.toWdlString})", ex)
        case (_, None) =>
          throw new VariableLookupException(s"Could not find a shard for scatter block with expression (${scatter.collection.toWdlString})")
      }
    }

    def handleDeclarationEvaluation(declaration: DeclarationInterface): Option[WdlValue] = {
      for {
        expression <- declaration.expression
        parentLookup = declaration.parent.map(_.lookupFunction(knownInputs, wdlFunctions, outputResolver, shards)).getOrElse(NoLookup)
        value = expression.evaluate(parentLookup, wdlFunctions).get
      } yield value
    }

    def handleCallEvaluation(call: Call): Option[WdlValue] = {
      this match {
          // Only use the shard number if the call is inside the scatter
        case s: Scatter if children.contains(call) => shards.get(s) map { shard =>
          outputResolver(call, Option(shard)) getOrElse {
            throw new VariableLookupException(s"Could not find outputs for call ${call.fullyQualifiedName} at shard $shard")
          }
        } orElse {
          throw new VariableLookupException(s"Could not find a shard for scatter block with expression (${s.collection.toWdlString})")
        }
        case _ => outputResolver(call, None).toOption
      }
    }

    def lookup(name: String): WdlValue = {
      val scopeResolvedValue = resolveVariable(name, relativeTo) flatMap {
        // First check if the variable has been provided in the known values
        case scope if knownInputs.contains(scope.fullyQualifiedName) => knownInputs.get(scope.fullyQualifiedName)
        case call: Call => handleCallEvaluation(call)
        case scatter: Scatter => handleScatterResolution(scatter)
        case declaration: DeclarationInterface if declaration.expression.isDefined => handleDeclarationEvaluation(declaration)
        case _ => None
      }

      if (scopeResolvedValue.isDefined) scopeResolvedValue.get
      else throw new VariableNotFoundException(name)
    }

    lookup
  }
}
