package wdl4s.wdl

import cats.data.Validated.{Invalid, Valid}
import lenthall.util.TryUtil
import wdl4s.parser.WdlParser._
import wdl4s.wdl.AstTools._
import wdl4s.wdl.exception.ValidationException
import wdl4s.wdl.expression.WdlFunctions
import wdl4s.wdl.types.WdlType
import wdl4s.wdl.values.WdlValue
import wdl4s.wom.callable.WorkflowDefinition
import wdl4s.wom.graph.{Graph, GraphNode}

import scala.language.postfixOps
import scala.util.Try

object WdlWorkflow {
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): WdlWorkflow = {
    if (ast.getName != AstNodeName.Workflow) {
      throw new UnsupportedOperationException(s"Expecting Workflow AST, got a ${ast.getName} AST")
    }
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val callNames = ast.findAsts(AstNodeName.Call).map {call =>
      ("Call", Option(call.getAttribute("alias")).getOrElse(call.getAttribute("task")))
    }
    val declarationNames = ast.findAsts(AstNodeName.Declaration).map {decl =>
      ("Declaration", decl.getAttribute("name"))
    }

    val workflowOutputsWildcards = ast.findAsts(AstNodeName.WorkflowOutputWildcard) map { wfOutput =>
      val wildcard = Option(wfOutput.getAttribute("wildcard")).map(_.sourceString).getOrElse("").nonEmpty
      val outputFqn = name + "." + wfOutput.getAttribute("fqn").sourceString
      WorkflowOutputWildcard(outputFqn, wildcard, wfOutput)
    }

    (callNames ++ declarationNames) groupBy { _._2.sourceString } foreach {
      case (_, terminals) if terminals.size > 1 =>
        val castedToTerminal = terminals map { case (astType, terminalAst) => (astType, terminalAst.asInstanceOf[Terminal])

        }

        throw new SyntaxError(wdlSyntaxErrorFormatter.multipleCallsAndHaveSameName(castedToTerminal))
      case _ =>
    }

    val meta = AstTools.wdlSectionToStringMap(ast, AstNodeName.Meta, wdlSyntaxErrorFormatter)
    val parameterMeta = AstTools.wdlSectionToStringMap(ast, AstNodeName.ParameterMeta, wdlSyntaxErrorFormatter)

    new WdlWorkflow(name, workflowOutputsWildcards, wdlSyntaxErrorFormatter, meta, parameterMeta, ast)
  }

  def buildWomGraph(wdlWorkflow: WdlWorkflow): Graph = {
    val graphNodes = wdlWorkflow.calls.foldLeft(Set.empty[GraphNode])({
      case (currentNodes, call) => currentNodes ++ call.womGraphInputNodes + call.womCallNode
    })

    Graph.validateAndConstruct(graphNodes) match {
      case Valid(wg) => wg.withDefaultOutputs
      case Invalid(errors) => throw ValidationException("Unable to validate graph", errors.map(new Exception(_)).toList)
    }
  }

  /**
    * Convert this WdlWorkflow into a wom.components.Workflow
    */
  def womWorkflowDefinition(wdlWorkflow: WdlWorkflow): WorkflowDefinition = {
    WorkflowDefinition(
      wdlWorkflow.unqualifiedName,
      buildWomGraph(wdlWorkflow),
      wdlWorkflow.meta,
      wdlWorkflow.parameterMeta,
      List.empty)
  }
}

case class WdlWorkflow(unqualifiedName: String,
                       workflowOutputWildcards: Seq[WorkflowOutputWildcard],
                       wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter,
                       meta: Map[String, String],
                       parameterMeta: Map[String, String],
                       ast: Ast) extends WdlCallable {

  /**
    * Convert this WdlWorkflow into a wom.components.Workflow
    */
  override lazy val womDefinition: WorkflowDefinition = WdlWorkflow.womWorkflowDefinition(this)

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
  override def resolveVariable(name: String, relativeTo: Scope = this): Option[WdlGraphNode] = {
    findCallByName(name) orElse findDeclarationByName(name) orElse findWorkflowOutputByName(name, relativeTo) orElse super.resolveVariable(name, relativeTo)
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
  def findCallByName(name: String): Option[WdlCall] = calls.find(_.unqualifiedName == name)

  /**
    * Declarations within the workflow scope (including inside scatters and ifs)
    */
  lazy val transitiveDeclarations = {
    def isValid(d: Declaration) = {
      d.parent match {
        case Some(w: WdlWorkflow) if w == this => true
        case Some(_: Scatter) => true
        case Some(_: If) => true
        case _ => false
      }
    }

    descendants collect {
      case declaration: Declaration if isValid(declaration) => declaration
    }
  }

  def findDeclarationByName(name: String): Option[Declaration] = {
    transitiveDeclarations.find(_.unqualifiedName == name)
  }

  def findWorkflowOutputByName(name: String, relativeTo: Scope) = {
    val leftOutputs = if (outputs.contains(relativeTo))
      outputs.dropRight(outputs.size - outputs.indexOf(relativeTo))
    else outputs
    leftOutputs.find(_.unqualifiedName == name)
  }

  lazy val hasEmptyOutputSection = workflowOutputWildcards.isEmpty && children.collect({ case o: WorkflowOutput => o }).isEmpty

  /**
   * All outputs for this workflow and their associated types
   *
   * @return a Map[FullyQualifiedName, WdlType] representing the union
   *         of all outputs from all `call`s within this workflow
   */
  lazy val expandedWildcardOutputs: Seq[WorkflowOutput] = {

    def toWorkflowOutput(output: DeclarationInterface, wdlType: WdlType) = {
      val locallyQualifiedName = output.parent map { parent => output.locallyQualifiedName(parent) } getOrElse {
        throw new RuntimeException(s"output ${output.fullyQualifiedName} has no parent Scope")
      }

      new WorkflowOutput(locallyQualifiedName, wdlType, WdlExpression.fromString(locallyQualifiedName), output.ast, Option(this))
    }

    def toWorkflowOutputs(scope: Scope) = {
      // Find out the number of parent scatters
      val outputs = scope match {
        case call: WdlCall => call.outputs
        case outputDeclaration: Output => Seq(outputDeclaration)
          // For non output declaration, don't return an array but return the raw value
        case otherDeclaration: DeclarationInterface => Seq(otherDeclaration)
        case _ => Seq.empty
      }

      outputs map { output => toWorkflowOutput(output, output.relativeWdlType(this)) }
    }

    // No outputs means all outputs
    val effectiveOutputWildcards = if (hasEmptyOutputSection) {
      calls map { call => WorkflowOutputWildcard(unqualifiedName + "." + call.unqualifiedName, wildcard = true, call.ast) } toSeq
    } else workflowOutputWildcards

    effectiveOutputWildcards flatMap { output =>
      // Prepend the namespace name in case it's an imported subworkflow so it can be resolved
      val outputFqn = namespace.unqualifiedName match {
        case empty if empty.isEmpty => output.fqn
        case alias => alias + "." + output.fqn
      }
      namespace.resolveCallOrOutputOrDeclaration(outputFqn) match {
        case Some(call: WdlCall) if output.wildcard && calls.contains(call) => toWorkflowOutputs(call)
        case Some(declaration: DeclarationInterface) if descendants.contains(declaration) => toWorkflowOutputs(declaration)
        case _ => throw new SyntaxError(wdlSyntaxErrorFormatter.badOldStyleWorkflowOutput(output.ast))
      }
    }
  }

  override lazy val outputs: Seq[WorkflowOutput] = expandedWildcardOutputs ++ children collect { case o: WorkflowOutput => o }

  override def toString = s"[Workflow $fullyQualifiedName]"

  def evaluateOutputs(knownInputs: WorkflowCoercedInputs,
                      wdlFunctions: WdlFunctions[WdlValue],
                      outputResolver: OutputResolver = NoOutputResolver,
                      shards: Map[Scatter, Int] = Map.empty[Scatter, Int]): Try[Map[WorkflowOutput, WdlValue]] = {

    val evaluatedOutputs = outputs.foldLeft(Map.empty[WorkflowOutput, Try[WdlValue]])((outputMap, output) => {
      val currentOutputs = outputMap collect {
        case (outputName, value) if value.isSuccess => outputName.fullyQualifiedName -> value.get
      }
      def knownValues = currentOutputs ++ knownInputs
      val lookup = lookupFunction(knownValues, wdlFunctions, outputResolver, shards, output)
      val coerced = output.requiredExpression.evaluate(lookup, wdlFunctions) flatMap output.wdlType.coerceRawValue
      val workflowOutput = output -> coerced

      outputMap + workflowOutput
    }) map { case (k, v) => k -> v }

    TryUtil.sequenceMap(evaluatedOutputs, "Failed to evaluate workflow outputs.\n")
  }
}

case class ReportableSymbol(fullyQualifiedName: FullyQualifiedName, wdlType: WdlType)
