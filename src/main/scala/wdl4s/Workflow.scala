package wdl4s

import lenthall.util.TryUtil
import wdl4s.AstTools.{AstNodeName, EnhancedAstNode}
import wdl4s.expression.WdlFunctions
import wdl4s.parser.WdlParser.{Ast, SyntaxError, Terminal}
import wdl4s.types.WdlType
import wdl4s.values.WdlValue

import scala.language.postfixOps
import scala.util.Try

object Workflow {
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Workflow = {
    if (ast.getName != AstNodeName.Workflow) {
      throw new UnsupportedOperationException(s"Expecting Workflow AST, got a ${ast.getName} AST")
    }
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val callNames = ast.findAsts(AstNodeName.Call).map {call =>
      Option(call.getAttribute("alias")).getOrElse(call.getAttribute("task"))
    }
    val workflowOutputsWildcards = ast.findAsts(AstNodeName.WorkflowOutputWildcard) map { wfOutput =>
      val wildcard = Option(wfOutput.getAttribute("wildcard")).map(_.sourceString).getOrElse("").nonEmpty
      val outputFqn = name + "." + wfOutput.getAttribute("fqn").sourceString
      WorkflowOutputWildcard(outputFqn, wildcard, wfOutput)
    }

    callNames groupBy { _.sourceString } foreach {
      case (_, terminals) if terminals.size > 1 =>
        throw new SyntaxError(wdlSyntaxErrorFormatter.multipleCallsAndHaveSameName(terminals.asInstanceOf[Seq[Terminal]]))
      case _ =>
    }

    val meta = AstTools.wdlSectionToStringMap(ast, AstNodeName.Meta, wdlSyntaxErrorFormatter)
    val parameterMeta = AstTools.wdlSectionToStringMap(ast, AstNodeName.ParameterMeta, wdlSyntaxErrorFormatter)

    new Workflow(name, workflowOutputsWildcards, wdlSyntaxErrorFormatter, meta, parameterMeta, ast)
  }
}

case class Workflow(unqualifiedName: String,
                    workflowOutputWildcards: Seq[WorkflowOutputWildcard],
                    wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter,
                    meta: Map[String, String],
                    parameterMeta: Map[String, String],
                    ast: Ast) extends Callable {

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
  override def resolveVariable(name: String, relativeTo: Scope = this): Option[GraphNode] = {
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
  def findCallByName(name: String): Option[Call] = calls.find(_.unqualifiedName == name)

  /**
    * Declarations within the workflow scope (including inside scatters and ifs)
    */
  lazy val transitiveDeclarations = {
    def isValid(d: Declaration) = {
      d.parent match {
        case Some(w: Workflow) if w == this => true
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
      
      new WorkflowOutput(output.locallyQualifiedName(this), wdlType, WdlExpression.fromString(locallyQualifiedName), output.ast, Option(this))
    }

    def toWorkflowOutputs(scope: Scope) = {
      // Find out the number of parent scatters
      val outputs = scope match {
        case call: Call => call.outputs
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
        case Some(call: Call) if output.wildcard && calls.contains(call) => toWorkflowOutputs(call)
        case Some(declaration: DeclarationInterface) if descendants.contains(declaration) => toWorkflowOutputs(declaration)
        case e => throw new SyntaxError(wdlSyntaxErrorFormatter.memberAccessReferencesBadCallOutput(output.ast))
      }
    }
  }
  
  override lazy val outputs: Seq[WorkflowOutput] = expandedWildcardOutputs ++ children collect { case o: WorkflowOutput => o }
  
  override def toString = s"[Workflow $fullyQualifiedName]"

  def evaluateOutputs(knownInputs: WorkflowCoercedInputs,
                      wdlFunctions: WdlFunctions[WdlValue],
                      outputResolver: OutputResolver = NoOutputResolver,
                      shards: Map[Scatter, Int] = Map.empty[Scatter, Int]): Try[Map[LocallyQualifiedName, WdlValue]] = {
    
    val evaluatedOutputs = outputs.foldLeft(Map.empty[WorkflowOutput, Try[WdlValue]])((outputMap, output) => {
      val currentOutputs = outputMap collect {
        case (outputName, value) if value.isSuccess => outputName.fullyQualifiedName -> value.get
      }
      def knownValues = currentOutputs ++ knownInputs
      val lookup = lookupFunction(knownValues, wdlFunctions, outputResolver, shards, output)
      val coerced = output.requiredExpression.evaluate(lookup, wdlFunctions) flatMap output.wdlType.coerceRawValue
      val workflowOutput = output -> coerced

      outputMap + workflowOutput
    }) map { case (k, v) => k.unqualifiedName -> v }

    TryUtil.sequenceMap(evaluatedOutputs, "Failed to evaluate workflow outputs.\n")
  }
}

case class ReportableSymbol(fullyQualifiedName: FullyQualifiedName, wdlType: WdlType)