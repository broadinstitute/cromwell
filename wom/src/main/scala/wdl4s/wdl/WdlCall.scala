package wdl4s.wdl

import cats.data.Validated.Valid
import cats.syntax.traverse._
import cats.instances.list._
import cats.syntax.cartesian._

import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.parser.WdlParser.{Ast, SyntaxError, Terminal}
import wdl4s.wdl.AstTools.EnhancedAstNode
import wdl4s.wdl.exception.{ValidationException, VariableLookupException, VariableNotFoundException}
import wdl4s.wdl.expression.WdlFunctions
import wdl4s.wdl.types.WdlOptionalType
import wdl4s.wdl.values.{WdlOptionalValue, WdlValue}
import wdl4s.wom.graph.CallNode.CallWithInputs
import wdl4s.wom.graph.{CallNode, GraphInputNode, GraphNodePort}
import cats.syntax.validated._
import lenthall.validation.ErrorOr.ShortCircuitingFlatMap

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

  private def buildWomNodeAndInputs(wdlCall: WdlCall): ErrorOr[CallWithInputs] = {

    val inputToOutputPortValidations: Iterable[ErrorOr[(String, GraphNodePort.OutputPort)]] = for {
      (inputName, expr) <- wdlCall.inputMappings
      variable <- expr.variableReferences
      parent <- wdlCall.parent
      node <- parent.resolveVariable(variable.terminal.sourceString)
      outputPort = outputPortFromNode(node, variable.terminalSubIdentifier)
      inputNameToOutputPortMapping = outputPort map { op => inputName -> op }
    } yield inputNameToOutputPortMapping

    val inputToOutputPortValidation: ErrorOr[List[(String, GraphNodePort.OutputPort)]] = inputToOutputPortValidations.toList.sequence

    // TODO: When WomExpressions wrap WdlExpressions, we can switch this to use expressionBasedInputs instead of portBasedInputs:
    val combinedValidation = (inputToOutputPortValidation |@| wdlCall.callable.womDefinition) map { (_, _) }

    for {
      combined <- combinedValidation
      (inputToOutputPort, callable) = combined
      callWithInputs <- CallNode.callWithInputs(wdlCall.alias.getOrElse(wdlCall.callable.unqualifiedName), callable, inputToOutputPort.toMap, Set.empty)
    } yield callWithInputs
  }

  //
  private def outputPortFromNode(node: WdlGraphNode, terminal: Option[Terminal]): ErrorOr[GraphNodePort.OutputPort] = {

    def findNamedOutputPort(name: String, graphOutputPorts: Set[GraphNodePort.OutputPort], terminalName: String): ErrorOr[GraphNodePort.OutputPort] = {
      graphOutputPorts.find(_.name == name) match {
        case Some(gop) => Valid(gop)
        case None => s"Cannot find an output port $name in $terminalName".invalidNel
      }
    }

    (node, terminal) match {
      case (wdlCall: WdlCall, Some(subTerminal)) =>
        for {
          graphOutputPorts <- wdlCall.womGraphOutputPorts
          namedPort <- findNamedOutputPort(subTerminal.sourceString, graphOutputPorts, "call " + wdlCall.unqualifiedName) // = graphOutputPorts.find(_.name == subTerminal.sourceString)
        } yield namedPort

        // TODO implement when declarations are in WOM
      case (_: Declaration, None) => "Declaration not yet supported in WOM".invalidNel
      case _ => s"Unsupported node $node and terminal $terminal".invalidNel
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

  private lazy val womCallNodeWithInputs = WdlCall.buildWomNodeAndInputs(this)

  lazy val womCallNode: ErrorOr[CallNode] = womCallNodeWithInputs.map(_.call)
  lazy val womGraphInputNodes: ErrorOr[Set[GraphInputNode]] = womCallNodeWithInputs.map(_.inputs)
  lazy val womGraphOutputPorts: ErrorOr[Set[GraphNodePort.OutputPort]] = womCallNode.map(_.outputPorts)

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