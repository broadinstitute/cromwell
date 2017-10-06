package wdl

import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.foldable._
import shapeless.Coproduct
import wdl4s.parser.WdlParser.{Ast, SyntaxError, Terminal}
import wdl.AstTools.EnhancedAstNode
import wdl.exception.{ValidationException, VariableLookupException, VariableNotFoundException}
import wdl.expression.WdlFunctions
import wdl.types.WdlOptionalType
import wdl.values._
import wom.callable.Callable
import wom.callable.Callable._
import wom.graph.CallNode._
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object WdlCall {
  def apply(ast: Ast,
            namespaces: Seq[WdlNamespace],
            tasks: Seq[WdlTask],
            workflows: Seq[WdlWorkflow],
            wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): WdlCall = {
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
      case task: WdlTask => WdlTaskCall(alias, task, callInputSectionMappings, ast)
      case workflow: WdlWorkflow => WdlWorkflowCall(alias, workflow, callInputSectionMappings, ast)
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

  private[wdl] def buildWomNodeAndInputs(wdlCall: WdlCall, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, GraphNodePort.OutputPort]) = {
    import lenthall.validation.ErrorOr._

    val callNodeBuilder = new CallNode.CallNodeBuilder()

    /*
      * Each input mapping gets its own ExpressionNode:
      * 
      * call my_task { input:
      *   input1 = "hi!"                -> ExpressionNode with no input port
      *   input3 = other_task.out + 2   -> ExpressionNode with an input port pointing to the output port of other_task.out
      * }
     */
    def expressionNodeMappings: ErrorOr[Map[LocalName, ExpressionNode]] = wdlCall.inputMappings traverse {
      case (inputName, wdlExpression) =>
        val identifier = wdlCall.womIdentifier.combine(inputName)
        WdlWomExpression.toExpressionNode(identifier, WdlWomExpression(wdlExpression, None), localLookup, outerLookup) map {
          LocalName(inputName) -> _
        }
    }

    /*
      * Fold over the input definitions and
      * 1) assign each input definition its InputDefinitionPointer
      * 2) if necessary, create a graph input node and assign its output port to the input definition
      * 
      * The InputDefinitionFold accumulates the input definition mappings, the create graph input nodes, and the expression nodes.
     */
    def foldInputDefinitions(expressionNodes: Map[LocalName, ExpressionNode], callable: Callable): InputDefinitionFold = {
      // Updates the fold with a new graph input node. Happens when an optional or required undefined input without an
      // expression node mapping is found
      def withGraphInputNode(inputDefinition: InputDefinition, graphInputNode: GraphInputNode) = {
        InputDefinitionFold(
          mappings = Map(inputDefinition -> Coproduct[InputDefinitionPointer](graphInputNode.singleOutputPort: OutputPort)),
          callInputPorts = Set(callNodeBuilder.makeInputPort(inputDefinition, graphInputNode.singleOutputPort)),
          newGraphInputNodes = Set(graphInputNode)
        )
      }

      callable.inputs.foldMap {
        // If there is an input mapping for this input definition, use that
        case inputDefinition if expressionNodes.contains(inputDefinition.localName) =>
          val expressionNode = expressionNodes(inputDefinition.localName)
          InputDefinitionFold(
            mappings = Map(inputDefinition -> expressionNode.inputDefinitionPointer),
            callInputPorts = Set(callNodeBuilder.makeInputPort(inputDefinition, expressionNode.singleExpressionOutputPort)),
            newExpressionNodes = Set(expressionNode)
          )

        // No input mapping, use the default expression
        case withDefault @ InputDefinitionWithDefault(_, _, expression) =>
          InputDefinitionFold(
            mappings = Map(withDefault -> Coproduct[InputDefinitionPointer](expression))
          )

        // No input mapping, required and we don't have a default value, create a new RequiredGraphInputNode
        // so that it can be satisfied via workflow inputs
        case required @ RequiredInputDefinition(n, womType) =>
          val identifier = wdlCall.womIdentifier.combine(n)
          withGraphInputNode(required, RequiredGraphInputNode(identifier, womType))

        // No input mapping, no default value but optional, create a OptionalGraphInputNode
        // so that it can be satisfied via workflow inputs
        case optional @ OptionalInputDefinition(n, womType) =>
          val identifier = wdlCall.womIdentifier.combine(n)
          withGraphInputNode(optional, OptionalGraphInputNode(identifier, womType))
      }
    }

    (expressionNodeMappings, wdlCall.callable.womDefinition) mapN {
      case (mappings, callable) =>
        callNodeBuilder.build(wdlCall.womIdentifier, callable, foldInputDefinitions(mappings, callable))
    }
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
sealed abstract class WdlCall(val alias: Option[String],
                              val callable: WdlCallable,
                              val inputMappings: Map[String, WdlExpression],
                              val ast: Ast) extends WdlGraphNodeWithInputs with WorkflowScoped {
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

    val declarationAttempts = callable.declarations.foldLeft[EvaluatedDeclarations](Map.empty)(doDeclaration)

    val (success, errors) = declarationAttempts partition {
      case (_, Success(_)) => true
      case _ => false
    }

    if (errors.nonEmpty) {
      val throwables = errors.toList map { _._2.failed.get }
      Failure(ValidationException(s"Input evaluation for Call $fullyQualifiedName failed.", throwables))
    } else {
      Success(success map { case (d, v) => d -> v.get })
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

      resolutions collectFirst { case Success(value) => value } getOrElse {
        resolutions.toList.flatMap({
          case Failure(_: VariableNotFoundException) => None
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

case class WdlTaskCall(override val alias: Option[String], task: WdlTask, override val inputMappings: Map[String, WdlExpression], override val ast: Ast) extends WdlCall(alias, task, inputMappings, ast) {
  override val callType = "call"
}
case class WdlWorkflowCall(override val alias: Option[String], calledWorkflow: WdlWorkflow, override val inputMappings: Map[String, WdlExpression], override val ast: Ast) extends WdlCall(alias, calledWorkflow, inputMappings, ast) {
  override val callType = "workflow"
}
