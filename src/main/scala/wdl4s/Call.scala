package wdl4s

import cats.data.NonEmptyList
import wdl4s.AstTools.EnhancedAstNode
import wdl4s.expression.WdlFunctions
import wdl4s.parser.WdlParser.{Ast, SyntaxError, Terminal}
import wdl4s.types.WdlOptionalType
import wdl4s.values.WdlValue

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object Call {
  def apply(ast: Ast,
            namespaces: Seq[WdlNamespace],
            tasks: Seq[Task],
            workflows: Seq[Workflow],
            wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Call = {
    val alias: Option[String] = ast.getAttribute("alias") match {
      case x: Terminal => Option(x.getSourceString)
      case _ => None
    }

    val taskName = ast.getAttribute("task").sourceString

    val callable = WdlNamespace.findCallable(taskName, namespaces, tasks ++ workflows) getOrElse {
      throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskName(ast, taskName))
    }

    val callInputSectionMappings = processCallInput(ast, wdlSyntaxErrorFormatter)

    callable match {
      case task: Task => new TaskCall(alias, task, callInputSectionMappings, ast)
      case workflow: Workflow => new WorkflowCall(alias, workflow, callInputSectionMappings, ast)
    }
  }

  private def processCallInput(ast: Ast,
                               wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Map[String, WdlExpression] = {
    AstTools.callInputSectionIOMappings(ast, wdlSyntaxErrorFormatter) map { a =>
      val key = a.getAttribute("key").sourceString
      val expression = new WdlExpression(a.getAttribute("value"))
      (key, expression)
    } toMap
  }
}

/**
 * Represents a `call` block in a WDL workflow.  Calls wrap tasks
 * and optionally provide a subset of the inputs for the task (as inputMappings member).
 * All remaining inputs to the task that are not declared in the `input` section
 * of the `call` block are called unsatisfiedInputs
 *
 * @param alias The alias for this call.  If two calls in a workflow use the same task
 *              then one of them needs to be aliased to a different name
 * @param callable The callable that this `call` will invoke
 * @param inputMappings A map of task-local input names and corresponding expression for the
 *                      value of those inputs
 */
sealed abstract class Call(val alias: Option[String],
                    val callable: Callable,
                    val inputMappings: Map[String, WdlExpression],
                    val ast: Ast) extends GraphNode with WorkflowScoped {
  val unqualifiedName: String = alias getOrElse callable.unqualifiedName
  
  def callType: String
  
  def toCallOutput(output: Output) = output match {
    case taskOutput: TaskOutput => CallOutput(this, taskOutput.copy(parent = Option(this)))
    case workflowOutput: WorkflowOutput => CallOutput(this, workflowOutput.copy(parent = Option(this)))
    case error => throw new Exception(s"Invalid output type ${error.getClass.getSimpleName}")
  }
  
  lazy val outputs: Seq[CallOutput] = callable.outputs map toCallOutput
  
  override def children: Seq[Scope] = super.children ++ outputs
  
  lazy val upstream: Set[GraphNode] = {
    val dependentNodes = for {
      expr <- inputMappings.values
      variable <- expr.variableReferences
      node <- parent.flatMap(_.resolveVariable(variable.sourceString))
    } yield node

    val firstScatterOrIf = ancestry.collectFirst({
      case s: Scatter with GraphNode => s
      case i: If with GraphNode => i
    })

    (dependentNodes ++ firstScatterOrIf.toSeq).toSet
  }

  lazy val downstream: Set[GraphNode] = {
    def expressions(node: GraphNode): Iterable[WdlExpression] = node match {
      case scatter: Scatter => Set(scatter.collection)
      case call: TaskCall => call.inputMappings.values
      case ifStatement: If => Set(ifStatement.condition)
      case declaration: Declaration => declaration.expression.toSet
      case _ => Set.empty
    }

    for {
      node <- namespace.descendants.collect({ 
        case n: GraphNode if n.fullyQualifiedNameWithIndexScopes != fullyQualifiedNameWithIndexScopes => n 
      })
      expression <- expressions(node)
      variable <- expression.variableReferences
      referencedNode = resolveVariable(variable.sourceString)
      if referencedNode == Option(this)
    } yield node
  }

  /**
   * Returns a Seq[WorkflowInput] representing the inputs to the call that are
   * needed before its command can be constructed. This excludes inputs that
   * are satisfied via the 'input' section of the Call definition.
   */
  def unsatisfiedInputs: Seq[WorkflowInput] = for {
    i <- declarations if !inputMappings.contains(i.unqualifiedName) && i.expression.isEmpty
  } yield WorkflowInput(i.fullyQualifiedName, i.wdlType)

  override def toString: String = s"[Call $fullyQualifiedName]"

  /**
    * The call is responsible for evaluating runtime inputs for its underlying task,
    * as the input value are provided for a specific call.
    * The returned value is a map from Declaration to WdlValue.
    * The keys int the return value are the task's declarations,
    * not the call's, as they will be used later for command instantiation
    * as well as output evaluation, which will both be performed by the task.
    */
  def evaluateTaskInputs(inputs: WorkflowCoercedInputs,
                         wdlFunctions: WdlFunctions[WdlValue],
                         outputResolver: OutputResolver = NoOutputResolver,
                         shards: Map[Scatter, Int] = Map.empty[Scatter, Int]): EvaluatedTaskInputs = {
    val declarationAttempts = callable.declarations map { declaration =>
      val lookup = lookupFunction(inputs, wdlFunctions, outputResolver, shards, relativeTo = declaration)
      val evaluatedDeclaration = Try(lookup(declaration.unqualifiedName))
      val coercedDeclaration = evaluatedDeclaration flatMap { declaration.wdlType.coerceRawValue(_) }
      
      declaration -> coercedDeclaration
    }
    
    val (success, errors) = declarationAttempts partition {
      case (_, Success(_)) => true
      case (d, Failure(_: VariableNotFoundException)) if d.wdlType.isInstanceOf[WdlOptionalType] => true
      case _ => false
    }
    
    if (errors.nonEmpty) {
      val errorsMessage = errors map { _._2.failed.get.getMessage }
      throw new ValidationException(
        s"Input evaluation for Call $fullyQualifiedName failed.", NonEmptyList.fromListUnsafe(errorsMessage.toList)
      )
    }

    val successfulDeclarations = success map { case (d, v) => v.toOption map { d -> _ } }
    successfulDeclarations.flatten.toMap
  }

  /**
    * Overrides the default lookup function to provide call specific resolution.
    */
  override def lookupFunction(inputs: WorkflowCoercedInputs,
                              wdlFunctions: WdlFunctions[WdlValue],
                              outputResolver: OutputResolver = NoOutputResolver,
                              shards: Map[Scatter, Int] = Map.empty[Scatter, Int],
                              relativeTo: Scope = this): String => WdlValue = {
    def lookup(name: String): WdlValue = {
      val inputMappingsWithMatchingName = Try(
        inputMappings.getOrElse(name, throw new Exception(s"Could not find $name in input section of call $fullyQualifiedName"))
      )

      val declarationsWithMatchingName = Try(
        declarations.find(_.unqualifiedName == name).getOrElse(throw new Exception(s"No declaration named $name for call $fullyQualifiedName"))
      )

      val inputMappingsLookup = for {
        inputExpr <- inputMappingsWithMatchingName
        parent <- Try(parent.getOrElse(throw new Exception(s"Call $unqualifiedName has no parent")))
        evaluatedExpr <- inputExpr.evaluate(parent.lookupFunction(inputs, wdlFunctions, outputResolver, shards, relativeTo), wdlFunctions)
      } yield evaluatedExpr

      val declarationLookup = for {
        declaration <- declarationsWithMatchingName
        inputsLookup <- Try(inputs.getOrElse(declaration.fullyQualifiedName, throw new Exception(s"No input for ${declaration.fullyQualifiedName}")))
      } yield inputsLookup

      val declarationExprLookup = for {
        declaration <- declarationsWithMatchingName
        declarationExpr <- Try(declaration.expression.getOrElse(throw new Exception(s"No expression defined for declaration ${declaration.fullyQualifiedName}")))
        evaluatedExpr <- declarationExpr.evaluate(lookupFunction(inputs, wdlFunctions, outputResolver, shards, relativeTo), wdlFunctions)
      } yield evaluatedExpr

      val taskParentResolution = for {
        parent <- Try(callable.parent.getOrElse(throw new Exception(s"Task ${callable.unqualifiedName} has no parent")))
        parentLookup <- Try(parent.lookupFunction(inputs, wdlFunctions, outputResolver, shards, relativeTo)(name))
      } yield parentLookup

      val resolutions = Seq(inputMappingsLookup, declarationLookup, declarationExprLookup, taskParentResolution)

      resolutions.collectFirst({ case s: Success[WdlValue] => s}) match {
        case Some(Success(value)) => value
        case None => throw new VariableNotFoundException(name)
      }
    }
    lookup
  }
}

case class TaskCall(override val alias: Option[String], task: Task, override val inputMappings: Map[String, WdlExpression], override val ast: Ast) extends Call(alias, task, inputMappings, ast) {
  override val callType = "call"
}
case class WorkflowCall(override val alias: Option[String], calledWorkflow: Workflow, override val inputMappings: Map[String, WdlExpression], override val ast: Ast) extends Call(alias, calledWorkflow, inputMappings, ast) {
  override val callType = "workflow"
}