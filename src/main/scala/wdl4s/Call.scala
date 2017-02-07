package wdl4s

import wdl4s.AstTools.EnhancedAstNode
import wdl4s.exception.{ValidationException, VariableLookupException, VariableNotFoundException}
import wdl4s.expression.WdlFunctions
import wdl4s.parser.WdlParser.{Ast, SyntaxError, Terminal}
import wdl4s.types.{WdlAnyType, WdlOptionalType}
import wdl4s.values.{WdlOptionalValue, WdlValue}

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
                    val ast: Ast) extends GraphNodeWithInputs with WorkflowScoped {
  val unqualifiedName: String = alias getOrElse callable.unqualifiedName
  
  def callType: String
  
  def toCallOutput(output: Output) = output match {
    case taskOutput: TaskOutput => CallOutput(this, taskOutput.copy(parent = Option(this)))
    case workflowOutput: WorkflowOutput => CallOutput(this, workflowOutput.copy(parent = Option(this)))
    case error => throw new Exception(s"Invalid output type ${error.getClass.getSimpleName}")
  }
  
  lazy val outputs: Seq[CallOutput] = callable.outputs map toCallOutput
  
  override def children: Seq[Scope] = super.children ++ outputs

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
                         shards: Map[Scatter, Int] = Map.empty[Scatter, Int]): Try[EvaluatedTaskInputs] = {

    type EvaluatedDeclarations = Map[Declaration, Try[WdlValue]]
    def doDeclaration(currentInputs: EvaluatedDeclarations, declaration: Declaration): EvaluatedDeclarations = {
      val newInputs = inputs ++ currentInputs.collect{
        case (decl, Success(value)) => decl.fullyQualifiedName -> value
      }
      val lookup = lookupFunction(newInputs, wdlFunctions, outputResolver, shards, relativeTo = declaration)
      val evaluatedDeclaration = Try(lookup(declaration.unqualifiedName))

      val coercedDeclaration: Try[WdlValue] = evaluatedDeclaration match {
        case Success(ed) => declaration.wdlType.coerceRawValue(ed)
        case Failure(_: VariableNotFoundException) if declaration.wdlType.isInstanceOf[WdlOptionalType] =>
          val innerType = declaration.wdlType.asInstanceOf[WdlOptionalType].memberType
          Success(WdlOptionalValue(innerType, None))
        case Failure(f) => Failure(f)
      }

      currentInputs + (declaration -> coercedDeclaration)
    }

    val declarationAttempts = callable.declarations.foldLeft[EvaluatedDeclarations](Map.empty)(doDeclaration _)

    val (success, errors) = declarationAttempts partition {
      case (_, Success(_)) => true
      case _ => false
    }

    if (errors.nonEmpty) {
      val throwables = errors.toList map { _._2.failed.get }
      Failure(ValidationException(s"Input evaluation for Call $fullyQualifiedName failed.", throwables))
    } else {
      Success(success map { case (d, v) => d -> v.get } toMap)
    }
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
        // Coerce the input into the declared type:
        declaration <- declarationsWithMatchingName
        coerced <- declaration.wdlType.coerceRawValue(evaluatedExpr)
      } yield coerced

      def unsuppliedDeclarationValue(declaration: Declaration) = declaration.wdlType match {
        case opt: WdlOptionalType => opt.none
        case _ => throw VariableNotFoundException(declaration)
      }

      val declarationLookup = for {
        declaration <- declarationsWithMatchingName
        inputsLookup <- Try(inputs.getOrElse(declaration.fullyQualifiedName, unsuppliedDeclarationValue(declaration)))
      } yield inputsLookup

      val declarationExprLookup = for {
        declaration <- declarationsWithMatchingName
        declarationExpr <- Try(declaration.expression.getOrElse(throw VariableNotFoundException(declaration)))
        evaluatedExpr <- declarationExpr.evaluate(lookupFunction(inputs, wdlFunctions, outputResolver, shards, relativeTo), wdlFunctions)
      } yield evaluatedExpr

      val taskParentResolution = for {
        parent <- Try(callable.parent.getOrElse(throw new Exception(s"Task ${callable.unqualifiedName} has no parent")))
        parentLookup <- Try(parent.lookupFunction(inputs, wdlFunctions, outputResolver, shards, relativeTo)(name))
      } yield parentLookup

      val resolutions = Seq(inputMappingsLookup, declarationExprLookup, declarationLookup, taskParentResolution)

      resolutions.collectFirst({ case Success(value) => value}) match {
        case Some(value) => value
        case None =>
          resolutions.toList.flatMap({
            case Failure(ex: VariableNotFoundException) => None
            case Failure(ex) => Option(ex) // Only take failures that are not VariableNotFoundExceptions
            case _ => None
          }) match {
            case Nil => throw VariableNotFoundException(name)
            case exs => throw new VariableLookupException(name, exs)
          }
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