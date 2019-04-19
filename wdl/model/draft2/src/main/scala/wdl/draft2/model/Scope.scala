package wdl.draft2.model

import common.Checked
import common.collections.EnhancedCollections._
import common.validation.Checked._
import wdl.draft2.model.exception.VariableNotFoundException
import wdl.draft2.model.expression.WdlFunctions
import wdl.draft2.model.values.WdlCallOutputsObject
import wdl.draft2.model.exception.{ScatterIndexNotFound, VariableLookupException}
import wdl.draft2.parser.WdlParser.Ast
import wom.graph.WomIdentifier
import wom.values.WomArray.WomArrayLike
import wom.values.WomValue

import scala.util.{Failure, Success, Try}
import scalax.collection.Graph
import scalax.collection.GraphEdge.DiEdge

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

  lazy val childGraphNodes: Set[WdlGraphNode] = children.toSet.filterByType[WdlGraphNode]
  lazy val childGraphNodesSorted: Checked[List[WdlGraphNode]] = {
    // We can't use a classic ordering because the upstream / downstream order is not total over the set of nodes
    // Instead we use topological sorting to guarantee that we process the nodes top to bottom
    val edges = for {
      child <- childGraphNodes
      upstreamNode <- child.upstreamAncestry if childGraphNodes.contains(upstreamNode)
    } yield DiEdge(upstreamNode, child)

    val graph = Graph.from[WdlGraphNode, DiEdge](childGraphNodes, edges)
    graph.topologicalSort match {
      case Left(_) =>
        val cycleNodes = childGraphNodes.filter(cgn => cgn.upstreamAncestry.contains(cgn)).map(_.toString).toList.sorted
        s"This workflow contains cyclic dependencies containing these edges: ${cycleNodes.mkString(", ")}".invalidNelCheck
      case Right(topologicalOrder) => topologicalOrder.toList.map(_.value).filter(childGraphNodes.contains).validNelCheck
    }
  }

  /**
    * Containing namespace
    */
  def namespace: WdlNamespace = _namespace
  private var _namespace: WdlNamespace = _
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
  private [wdl] def ancestrySafe: Seq[Scope] = parent match {
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
  lazy val calls: Set[WdlCall] = descendants.collect({ case c: WdlCall => c })

  lazy val taskCalls: Set[WdlTaskCall] = calls collect { case c: WdlTaskCall => c }

  lazy val workflowCalls: Set[WdlWorkflowCall] = calls collect { case c: WdlWorkflowCall => c }

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
  lazy val fullyQualifiedName = {
    (ancestry.filter(_.appearsInFqn).map(_.unqualifiedName).reverse :+ unqualifiedName).mkString(".")
  }
  
  lazy val womIdentifier = {
    /* Limit the fully qualified name to the parent workflow, if it exists.
     * The reason for this is if this scope comes from an imported namespace,
     * the name of the namespace will be part of the FQN. While useful to construct the hierarchy and dependencies,
     * this is not desirable when being used by underlying engines.
    */
    val womFullyQualifiedName = ancestry.collectFirst {
      case workflow: WdlWorkflow => locallyQualifiedName(workflow)
    } getOrElse fullyQualifiedName

    WomIdentifier(unqualifiedName, womFullyQualifiedName)
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
    val otherAncestry = other.ancestrySafe
    ancestrySafe find { otherAncestry.contains(_) }
  }

  def closestCommonScatter(other: Scope): Option[Scatter] = {
    val otherAncestry = other.ancestrySafe
    ancestrySafe collectFirst  {
      case s: Scatter if otherAncestry.contains(s) => s
    }
  }

  /**
    * Performs scope resolution starting from this scope and walking up the lexical hierarchy
    * until it finds a GraphNode with the `name` as its unqualifiedName
    */
  def resolveVariable(name: String, relativeTo: Scope = this, ignoreLocal: Boolean = false): Option[WdlGraphNode] = {
    val siblingScopes = if (children.contains(relativeTo))
    // For declarations, only resolve to declarations that are lexically before this declaration
      children.dropRight(children.size - children.indexOf(relativeTo) )
    else children

    val siblingCallOutputs = siblingScopes flatMap {
      case c: WdlTaskCall => c.outputs
      case _ => Seq.empty[CallOutput]
    }

    val localLookup =
      if (ignoreLocal)
        Seq.empty
      else (siblingScopes ++ siblingCallOutputs) collect {
        case d: Declaration if d.unqualifiedName == name => d
        case c: WdlTaskCall if c.unqualifiedName == name => c
        case co: CallOutput if co.unqualifiedName == name => co
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
    * @return String => WomValue lookup function rooted at `scope`
    * @throws VariableNotFoundException If no errors occurred, but also `name` didn't resolve to any value
    * @throws VariableLookupException if anything else goes wrong in looking up a value for `name`
    */
  def lookupFunction(knownInputs: WorkflowCoercedInputs,
                     wdlFunctions: WdlFunctions[WomValue],
                     outputResolver: OutputResolver = NoOutputResolver,
                     shards: Map[Scatter, Int] = Map.empty[Scatter, Int],
                     relativeTo: Scope = this): String => WomValue = {

    def handleScatterResolution(scatter: Scatter): Try[WomValue] = {
      // This case will happen if `name` references a Scatter.item (i.e. `x` in expression scatter(x in y) {...})
      val evaluatedCollection = scatter.collection.evaluate(scatter.lookupFunction(knownInputs, wdlFunctions, outputResolver, shards), wdlFunctions)
      val scatterShard = shards.get(scatter)

      (evaluatedCollection, scatterShard) match {
        case (Success(WomArrayLike(array)), Some(shard)) if 0 <= shard && shard < array.value.size =>
          array.value.lift(shard) match {
            case Some(v) => Success(v)
            case None => Failure(new VariableLookupException(s"Could not find value for shard index $shard in scatter collection $array"))
          }
        case (Success(WomArrayLike(array)), Some(shard)) =>
          Failure(new VariableLookupException(s"Scatter expression (${scatter.collection.toWomString}) evaluated to an array of ${array.value.size} elements, but element $shard was requested."))
        case (Success(_: WomArrayLike), None) =>
          Failure(ScatterIndexNotFound(s"Could not find the shard mapping to this scatter ${scatter.fullyQualifiedName}"))
        case (Success(value: WomValue), _) =>
          Failure(new VariableLookupException(s"Expected scatter expression (${scatter.collection.toWomString}) to evaluate to an Array.  Instead, got a $value"))
        case (failure @ Failure(_), _) => failure
        case (_, None) =>
          Failure(new VariableLookupException(s"Could not find a shard for scatter block with expression (${scatter.collection.toWomString})"))
      }
    }

    def fromOutputs(node: WdlGraphNode) = {
      def withShard(s: Scatter) = {
        shards.get(s) map { shard =>
          outputResolver(node, Option(shard))
        } getOrElse {
          Failure(ScatterIndexNotFound(s"Could not find a shard for scatter block with expression (${s.collection.toWomString})"))
        }
      }

      this match {
        case s: Scatter if descendants.contains(node) => withShard(s)
        case other =>
          other.closestCommonScatter(node) match {
            case Some(s: Scatter) => withShard(s)
            case _ => outputResolver(node, None)
          }
      }
    }

    def handleDeclarationEvaluation(declaration: DeclarationInterface): Try[WomValue] = {
      def declarationExcludingLookup(lookup: String => WomValue): String => WomValue = name => {
        if (name.equals(declaration.unqualifiedName)) throw new IllegalArgumentException(s"Declaration for '$name' refers to its own value")
        else lookup(name)
      }

      def evaluate =
        declaration.expression match {
          case Some(e) =>
            val parentLookup = declaration.parent.map(_.lookupFunction(knownInputs, wdlFunctions, outputResolver, shards)).getOrElse(NoLookup)
            e.evaluate(declarationExcludingLookup(parentLookup), wdlFunctions)
          case None =>
            Failure(new VariableLookupException(s"Declaration ${declaration.fullyQualifiedName} does not have an expression"))
        }

      fromOutputs(declaration) recoverWith { case _ => evaluate }
    }

    def handleCallEvaluation(call: WdlCall): Try[WomValue] = fromOutputs(call)

    def lookup(name: String): WomValue = {
      val scopeResolvedValue = resolveVariable(name, relativeTo) map {
        // First check if the variable has been provided in the known values
        case scope if knownInputs.contains(scope.fullyQualifiedName) =>
          knownInputs.get(scope.fullyQualifiedName) map Success.apply getOrElse {
            Failure(new VariableLookupException(s"Could not find value in inputs map."))
          }
        case callOutput: CallOutput => handleCallEvaluation(callOutput.call) flatMap {
          case outputs: WdlCallOutputsObject => outputs.outputs.get(callOutput.unqualifiedName).map(Success(_)).getOrElse(Failure(new Exception(s"No output ${callOutput.unqualifiedName} found in ${callOutput.call.unqualifiedName}'s outputs")))
          case other => Failure(new Exception(s"Call outputs unexpectedly evaluated to a ${other.womType.stableName}"))
        }
        case call: WdlCall => handleCallEvaluation(call)
        case scatter: Scatter => handleScatterResolution(scatter)
        case declaration: DeclarationInterface if declaration.expression.isDefined => handleDeclarationEvaluation(declaration)
        case scope => Failure(new VariableLookupException(s"Variable $name resolved to scope ${scope.fullyQualifiedName} but cannot be evaluated."))
      } getOrElse {
        Failure(VariableNotFoundException(name))
      }

      scopeResolvedValue.get
    }

    lookup
  }
}
