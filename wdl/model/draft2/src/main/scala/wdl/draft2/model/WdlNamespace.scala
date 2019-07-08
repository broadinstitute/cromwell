package wdl.draft2.model

import java.nio.file.{Path, Paths}

import better.files._
import common.util.TryUtil
import wdl.draft2.Draft2ResolvedImportBundle
import wdl.draft2.model.AstTools.{AstNodeName, EnhancedAstNode}
import wdl.draft2.model.command.ParameterCommandPart
import wdl.draft2.model.expression.WdlStandardLibraryFunctions
import wdl.draft2.model.types.WdlCallOutputsObjectType
import wdl.draft2.model.exception._
import wdl.draft2.model.expression.{NoFunctions, WdlStandardLibraryFunctionsType}
import wdl.draft2.model.types.WdlNamespaceType
import wdl.draft2.parser.WdlParser
import wdl.draft2.parser.WdlParser._
import wom.ResolvedImportRecord
import wom.core._
import wom.types._
import wom.values.WomValue

import scala.collection.JavaConverters._
import scala.collection.mutable
import scala.util.{Failure, Success, Try}


/**
  * Represents a parsed WDL file
  */
sealed trait WdlNamespace extends WomValue with Scope {
  final val womType = WdlNamespaceType
  def ast: Ast
  def resource = ast.findFirstTerminal.map(_.getResource).getOrElse("NONE")
  def importedAs: Option[String] // Used when imported with `as`
  def imports: Seq[Import]
  def namespaces: Seq[WdlNamespace]
  /** Return this `WdlNamespace` and its child `WdlNamespace`s recursively. */
  def allNamespacesRecursively: List[WdlNamespace] = {
    this +: (namespaces.toList flatMap { _.allNamespacesRecursively })
  }
  def tasks: Seq[WdlTask]
  def workflows: Seq[WdlWorkflow]
  def terminalMap: Map[Terminal, WorkflowSource]
  def findTask(name: String): Option[WdlTask] = tasks.find(_.name == name)
  override def unqualifiedName: LocallyQualifiedName = importedAs.getOrElse("")
  override def appearsInFqn: Boolean = importedAs.isDefined
  override def namespace: WdlNamespace = this
  def resolve(fqn: FullyQualifiedName): Option[Scope] = {
    (descendants + this).find(d => d.fullyQualifiedName == fqn || d.fullyQualifiedNameWithIndexScopes == fqn)
  }
  def resolveCallOrOutputOrDeclaration(fqn: FullyQualifiedName): Option[Scope] = {
    val callsAndOutputs = descendants collect {
      case c: WdlCall => c
      case d: Declaration => d
      case o: TaskOutput => o
      case o: CallOutput => o
    }
    callsAndOutputs.find(d => d.fullyQualifiedName == fqn || d.fullyQualifiedNameWithIndexScopes == fqn)
  }
  def sourceString: String
  def importUri: Option[String] // URI under which this namespace was imported
  def resolvedImportRecords: Set[ResolvedImportRecord]
}

/**
  * A WdlNamespace which doesn't have a locally defined Workflow.
  */
case class WdlNamespaceWithoutWorkflow(importedAs: Option[String],
                                       imports: Seq[Import],
                                       namespaces: Seq[WdlNamespace],
                                       tasks: Seq[WdlTask],
                                       terminalMap: Map[Terminal, WorkflowSource],
                                       ast: Ast,
                                       sourceString: String,
                                       importUri: Option[String] = None,
                                       resolvedImportRecords: Set[ResolvedImportRecord] = Set.empty[ResolvedImportRecord]) extends WdlNamespace {
  val workflows = Seq.empty[WdlWorkflow]

}

/**
  * A WdlNamespace which has exactly one workflow defined.
  */
case class WdlNamespaceWithWorkflow(importedAs: Option[String],
                                    workflow: WdlWorkflow,
                                    imports: Seq[Import],
                                    namespaces: Seq[WdlNamespace],
                                    tasks: Seq[WdlTask],
                                    terminalMap: Map[Terminal, WorkflowSource],
                                    wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter,
                                    ast: Ast,
                                    sourceString: String,
                                    importUri: Option[String] = None,
                                    resolvedImportRecords: Set[ResolvedImportRecord] = Set.empty[ResolvedImportRecord]) extends WdlNamespace {

  override val workflows = Seq(workflow)

  override def toString: String = s"[WdlNamespace importedAs=$importedAs]"

  /**
    * Some declarations need a value from the user and some have an expression attached to them.
    * For the declarations that have an expression attached to it already, evaluate the expression
    * and return the value. Only evaluates workflow level declarations. Other declarations will be evaluated at runtime.
    */
  def staticDeclarationsRecursive(userInputs: WorkflowCoercedInputs, wdlFunctions: WdlStandardLibraryFunctions): Try[WorkflowCoercedInputs] = {
    import common.exception.Aggregation._

    def evalDeclaration(accumulated: Map[FullyQualifiedName, Try[WomValue]], current: Declaration): Map[FullyQualifiedName, Try[WomValue]] = {
      current.expression match {
        case Some(expr) =>
          val successfulAccumulated = accumulated.collect({ case (k, v) if v.isSuccess => k -> v.get })
          val value = expr.evaluate(current.lookupFunction(successfulAccumulated ++ userInputs, wdlFunctions, NoOutputResolver, Map.empty[Scatter, Int]), wdlFunctions)
          val correctlyCoerced = value flatMap current.womType.coerceRawValue
          accumulated + (current.fullyQualifiedName -> correctlyCoerced)
        case None => accumulated
      }
    }

    def evalScope: Map[FullyQualifiedName, Try[WomValue]] = {
      val workflowDeclarations = children.collect({ case w: WdlWorkflow => w.declarations }).flatten
      (declarations ++ workflowDeclarations).foldLeft(Map.empty[FullyQualifiedName, Try[WomValue]])(evalDeclaration)
    }

    val filteredExceptions: Set[Class[_ <: Throwable]] = Set(classOf[OutputVariableLookupException], classOf[ScatterIndexNotFound])

    // Filter out declarations for which evaluation failed because a call output variable could not be resolved, or a shard could not be found,
    // as this method is meant for pre-execution validation
    val filtered = evalScope filterNot {
      case (_, Failure(ex)) => ex.flatten forall { ex => filteredExceptions.contains(ex.getClass) }
      case _ => false
    } map {
      case (name, Failure(f)) => name -> Failure(ValidationException(name, List(f)))
      case other => other
    }

    TryUtil.sequenceMap(filtered, "Could not evaluate workflow declarations")
  }
}

/**
  * Main interface into the `wdl4s` package.
  *
  * Example usage
  *
  * {{{
  * val namespace = WdlNamespace.loadUsingPath(Paths.get("/path/to/file.wdl"), None, None)
  * namespace.workflow.calls foreach { call =>
  *   println(call)
  * }
  * }}}
  */
object WdlNamespace {

  val WorkflowResourceString = "string"

  def loadUsingPath(wdlFile: Path, resource: Option[String], importResolver: Option[Seq[Draft2ImportResolver]]): Try[WdlNamespace] = {
    load(readFile(wdlFile), resource.getOrElse(wdlFile.toString), importResolver.getOrElse(Seq(fileResolver)), None)
  }

  def loadUsingSource(workflowSource: WorkflowSource, resource: Option[String], importResolver: Option[Seq[Draft2ImportResolver]]): Try[WdlNamespace] = {
    load(workflowSource, resource.getOrElse(WorkflowResourceString), importResolver.getOrElse(Seq(fileResolver)), None)
  }

  private def load(workflowSource: WorkflowSource,
                   resource: String,
                   importResolver: Seq[Draft2ImportResolver],
                   importedAs: Option[String],
                   resolvedImportRecords: Set[ResolvedImportRecord] = Set.empty[ResolvedImportRecord],
                   root: Boolean = true): Try[WdlNamespace] = Try {
    val maybeAst = Option(AstTools.getAst(workflowSource, resource))

    maybeAst match {
      case Some(ast) =>
        WdlNamespace(ast, resource, workflowSource, importResolver, importedAs, resolvedImportRecords, root = root)
      case None =>
        throw new IllegalArgumentException("Could not build AST from workflow source. Source is empty or contains only comments and whitespace.")
    }
  }

  def apply(ast: Ast,
            uri: String,
            source: WorkflowSource,
            importResolvers: Seq[Draft2ImportResolver],
            namespaceName: Option[String],
            resolvedImportRecords: Set[ResolvedImportRecord],
            root: Boolean = false): WdlNamespace = {

    val imports = for {
      importAst <- Option(ast).map(_.getAttribute("imports")).toSeq
      importStatement <- importAst.astListAsVector.map(Import(_))
    } yield importStatement

    def tryResolve(str: String, remainingResolvers: Seq[Draft2ImportResolver], errors: List[Throwable]): Draft2ResolvedImportBundle = {
      remainingResolvers match {
        case resolver :: tail =>
          Try(resolver(str)) match {
            case Success(s) => s
            case Failure(f) => tryResolve(str, tail, errors :+ f)
          }
        case Nil =>
          errors match {
            case Nil => throw new UnsatisfiedInputException("Failed to import workflow, no import sources provided.")
            case _ => throw ValidationException(s"Failed to import workflow $str.", errors)
          }
      }
    }

    // Translates all import statements to sub-namespaces
    val namespaces: Seq[WdlNamespace] = for {
      imp <- imports
      draft2ResolvedImportBundle = tryResolve(imp.uri, importResolvers, List.empty)
    } yield WdlNamespace.load(
      draft2ResolvedImportBundle.source,
      imp.uri,
      importResolvers,
      Option(imp.namespaceName),
      resolvedImportRecords + draft2ResolvedImportBundle.resolvedImportRecord,
      root = false
    ).get

    // Map of Terminal -> WDL Source Code so the syntax error formatter can show line numbers
    val terminalMap = AstTools.terminalMap(ast, source)
    val combinedTerminalMap = ((namespaces map { x => x.terminalMap }) ++ Seq(terminalMap)) reduce (_ ++ _)
    val wdlSyntaxErrorFormatter = WdlSyntaxErrorFormatter(combinedTerminalMap)

    val topLevelAsts = ast.getAttribute("body").astListAsVector.collect({ case a: Ast => a })

    // All `task` definitions of primary workflow.
    val topLevelTasks: Seq[WdlTask] = for {
      taskAst <- topLevelAsts if taskAst.getName == AstNodeName.Task
    } yield WdlTask(taskAst, wdlSyntaxErrorFormatter)

    val workflows: Seq[WdlWorkflow] = for {
      workflowAst <- topLevelAsts if workflowAst.getName == AstNodeName.Workflow
    } yield WdlWorkflow(workflowAst, wdlSyntaxErrorFormatter)

    // Build scope tree recursively
    val scopeIndexes: mutable.Map[Class[_ <: Scope], Int] = mutable.HashMap.empty.withDefaultValue(-1)

    def getScope(scopeAst: Ast, parent: Option[Scope]): Scope = {
      val scope = scopeAst.getName match {
        case AstNodeName.Call => WdlCall(scopeAst, namespaces, topLevelTasks, workflows, wdlSyntaxErrorFormatter)
        case AstNodeName.Workflow => WdlWorkflow(scopeAst, wdlSyntaxErrorFormatter)
        case AstNodeName.Declaration => Declaration(scopeAst, wdlSyntaxErrorFormatter, parent)
        case AstNodeName.Scatter =>
          scopeIndexes(classOf[Scatter]) += 1
          Scatter(scopeAst, scopeIndexes(classOf[Scatter]))
        case AstNodeName.If =>
          scopeIndexes(classOf[If]) += 1
          If(scopeAst, scopeIndexes(classOf[If]))
        case AstNodeName.Output => TaskOutput(scopeAst, wdlSyntaxErrorFormatter, parent)
        case AstNodeName.WorkflowOutputDeclaration => WorkflowOutput(scopeAst, wdlSyntaxErrorFormatter, parent)
      }

      scope.children = getChildren(scopeAst, Option(scope))
      scope.children.foreach(_.parent = scope)
      scope
    }

    def getChildren(scopeAst: Ast, scope: Option[Scope]): Seq[Scope] = {
      val ScopeAstNames = Seq(
        AstNodeName.Call, AstNodeName.Workflow, AstNodeName.Namespace,
        AstNodeName.Scatter, AstNodeName.If, AstNodeName.Declaration
      )

      def getScopeAsts(root: Ast, astAttribute: String): Seq[Ast] = {
        root.getAttribute(astAttribute).astListAsVector.collect({ case a: Ast if ScopeAstNames.contains(a.getName) => a })
      }

      def getTaskInputsOutputs(ast: Ast) = {
        val inputDeclarations = getScopeAsts(ast, "declarations").map(getScope(_, scope))
        val outputDeclarations = ast.findAsts(AstNodeName.Output).map(getScope(_, scope))
        inputDeclarations ++ outputDeclarations
      }

      def getWorkflowOutputs(ast: Ast) = ast.findAsts(AstNodeName.WorkflowOutputDeclaration).map(getScope(_, scope))

      scopeAst.getName match {
        case AstNodeName.Task => getTaskInputsOutputs(scopeAst)
        case AstNodeName.Declaration | AstNodeName.Output | AstNodeName.WorkflowOutputDeclaration => Seq.empty[Scope]
        case AstNodeName.Call =>
          val referencedTask = findCallable(scopeAst.getAttribute("task").sourceString, namespaces, topLevelTasks ++ workflows)
          referencedTask match {
            case Some(task: WdlTask) =>
              getScopeAsts(task.ast, "declarations").map(d => getScope(d, scope))
            case Some(workflow: WdlWorkflow) =>
              workflow.ast.getAttribute("body").astListAsVector collect {
                case declaration: Ast if declaration.getName == AstNodeName.Declaration => getScope(declaration, scope)
              }
            case _ => Seq.empty[Scope]
          }
        case AstNodeName.Scatter | AstNodeName.If | AstNodeName.Namespace =>
          getScopeAsts(scopeAst, "body").map(getScope(_, scope))
        case AstNodeName.Workflow =>
          getScopeAsts(scopeAst, "body").map(getScope(_, scope)) ++ getWorkflowOutputs(scopeAst)
      }
    }

    val topLevelDeclarationScopes = for {
      ast <- topLevelAsts
      if ast.getName != AstNodeName.Task && ast.getName != AstNodeName.Workflow
    } yield ast

    val children = topLevelTasks ++ namespaces ++ workflows ++ topLevelDeclarationScopes.map(ast => getScope(ast, parent = None))

    val newResolvedImportRecordsSet = resolvedImportRecords ++ namespaces.flatMap(_.resolvedImportRecords).toSet

    val namespace = workflows match {
      case Nil => WdlNamespaceWithoutWorkflow(
        importedAs = namespaceName,
        imports = imports,
        namespaces = namespaces,
        tasks = topLevelTasks,
        terminalMap = terminalMap,
        ast = ast,
        sourceString = source,
        importUri = Option(uri),
        resolvedImportRecords = newResolvedImportRecordsSet
      )
      case Seq(w) => WdlNamespaceWithWorkflow(
        ast = ast,
        workflow = w,
        namespace = namespaceName,
        imports = imports,
        namespaces = namespaces,
        tasks = topLevelTasks,
        terminalMap = terminalMap,
        wdlSyntaxErrorFormatter = wdlSyntaxErrorFormatter,
        sourceString = source,
        importUri = Option(uri),
        resolvedImportRecords = newResolvedImportRecordsSet
      )
      case _ => throw new SyntaxError(wdlSyntaxErrorFormatter.tooManyWorkflows(ast.findAsts(AstNodeName.Workflow).asJava))
    }

    // Write-once var setting for parent/child relationships
    def descendants(scope: Scope): Seq[Scope] = {
      val children = scope.children
      val childDescendants = scope.children.flatMap({
        case _: WdlNamespace => Seq.empty
        case s => descendants(s)
      })
      children ++ childDescendants
    }

    namespace.children = children
    namespace.children.foreach(_.parent = namespace)

    topLevelTasks foreach { task =>
      task.children = getChildren(task.ast, Option(task))
      task.children.foreach(_.parent = task)
    }

    workflows foreach { workflow =>
      workflow.children = getChildren(workflow.ast, Option(workflow))
      workflow.children.foreach(_.parent = workflow)
    }

    descendants(namespace).foreach(_.namespace = namespace)

    // SYNTAX CHECKS
    val callInputSectionErrors = namespace.descendants.collect({ case c: WdlCall => c }).flatMap(
      validateCallInputSection(_, wdlSyntaxErrorFormatter)
    )

    val workflowOutputErrors = workflows flatMap { _.workflowCalls map { _.calledWorkflow } } collect {
      case calledWorkflow if calledWorkflow.workflowOutputWildcards.nonEmpty =>
        new SyntaxError(
          s"""Workflow ${calledWorkflow.unqualifiedName} is used as a sub workflow but has outputs declared with a deprecated syntax not compatible with sub workflows.
             |To use this workflow as a sub workflow please update the workflow outputs section to the latest WDL specification.
             |See https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#outputs""".stripMargin
        )
    }

    /*
     * Used to filter task declaration nodes whose parents are subworkflow calls. Subworkflows will cause task declarations to appear twice among
     * the namespace's descendants: once with a parent that is the containing static workflow and again with a parent that is the dynamic workflow
     * call. Looking for variable references when the parent is a workflow call will fail if the variable references a call output since only
     * declarations are made children of workflow calls. This variable reference check is also redundant to the check that will happen with the
     * copy of the declaration that has a static workflow parent.
     */
    def parentIsNotAWorkflowCall(declaration: DeclarationInterface): Boolean = declaration.parent match {
      case Some(_: WdlWorkflowCall) => false
      case _ => true
    }

    val declarationErrors = for {
      descendant <- namespace.descendants
      declaration <- getDecls(descendant)
      error <- validateDeclaration(declaration, wdlSyntaxErrorFormatter) if parentIsNotAWorkflowCall(declaration)
    } yield error

    val runtimeErrors = for {
      task <- namespace.tasks
      runtime <- task.runtimeAttributes.attrs
      error <- validateRuntime(runtime._2, task, wdlSyntaxErrorFormatter)
    } yield error

    val scatterErrors = for {
      scatter <- namespace.descendants.collect { case sc: Scatter => sc }
      expression = scatter.collection
      badVariable <- referencesToAbsentValues(scatter, expression)
    } yield new SyntaxError(wdlSyntaxErrorFormatter.scatterCollectionContainsInvalidVariableReference(scatter, badVariable))

    def scopeNameAndTerminal(scope: Scope): (String, Terminal) = {
      scope match {
        case ns: WdlNamespace => ("Namespace", imports.find(_.uri == ns.resource).get.namespaceTerminal)
        case s: Scope => (s.getClass.getSimpleName, s.ast.findFirstTerminal.get)
      }
    }

    case class ScopeAccumulator(accumulated: Seq[Scope] = Seq.empty, errors: Seq[String] = Seq.empty)

    def lookForDuplicates(scopes: Traversable[Scope]) = {
      scopes.foldLeft(ScopeAccumulator()) { (acc, cur) =>
        val possibleError = acc.accumulated.find(_.unqualifiedName == cur.unqualifiedName) map { duplicate =>
          val (dupName, dupTerminal) = scopeNameAndTerminal(duplicate)
          val (curName, curTerminal) = scopeNameAndTerminal(cur)
          wdlSyntaxErrorFormatter.twoSiblingScopesHaveTheSameName(
            dupName, dupTerminal, curName, curTerminal
          )
        }
        ScopeAccumulator(acc.accumulated :+ cur, acc.errors ++ possibleError.toSeq)
      }
    }

    val scopeDuplicationErrors = (namespace.descendants + namespace) collect {
      case scope if scope.namespace == namespace => lookForDuplicates(scope.children)
    }

    val expandedWorkflowOutputsDuplicationErrors = {
      (namespace.descendants + namespace) collect { case workflow: WdlWorkflow => lookForDuplicates(workflow.outputs) }
    }

    val accumulatedErrors = expandedWorkflowOutputsDuplicationErrors ++ scopeDuplicationErrors

    val duplicateSiblingScopeNameErrors = accumulatedErrors.flatMap(_.errors).map(new SyntaxError(_)).toSeq

    val taskCommandReferenceErrors = for {
      task <- namespace.tasks
      param <- task.commandTemplate.collect({ case p: ParameterCommandPart => p })
      variable <- param.expression.variableReferences(task)
      if !task.declarations.map(_.unqualifiedName).contains(variable.terminal.getSourceString)
    } yield new SyntaxError(wdlSyntaxErrorFormatter.commandExpressionContainsInvalidVariableReference(task.ast.getAttribute("name").asInstanceOf[Terminal], variable.terminal))

    val all = workflowOutputErrors ++ declarationErrors ++ runtimeErrors ++ scatterErrors ++ callInputSectionErrors ++ taskCommandReferenceErrors ++ duplicateSiblingScopeNameErrors

    all.sortWith({ case (l, r) => l.getMessage < r.getMessage }) match {
      case s: Seq[SyntaxError] if s.nonEmpty => throw s.head
      case _ =>
    }

    namespace
  }

  private def getDecls(scope: Scope): Seq[DeclarationInterface] = {
    scope match {
      case t: WdlTask => t.declarations ++ t.outputs
      case w: WdlWorkflow => w.declarations ++ w.outputs
      case s => s.declarations
    }
  }

  private def declarationName(declarationAst: Ast): Terminal = declarationAst.getAttribute("name").asInstanceOf[Terminal]

  /**
    * Determine the list of references in this expression to values which were never declared
    */
  private def referencesToAbsentValues(container: Scope, expression: WdlExpression): Iterable[Terminal] =
    expression.variableReferences(container) collect { case variable if container.resolveVariable(variable.terminal.sourceString).isEmpty => variable.terminal }

  private def validateDeclaration(declaration: DeclarationInterface, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Seq[SyntaxError] = {
    val invalidVariableReferences = for {
      expr <- declaration.expression.toSeq
      variable <- referencesToAbsentValues(declaration, expr)
    } yield new SyntaxError(wdlSyntaxErrorFormatter.declarationContainsReferenceToAbsentValue(
      declaration.parent,
      variable
    ))

    val typeMismatches = typeCheckDeclaration(declaration, wdlSyntaxErrorFormatter).toSeq

    invalidVariableReferences ++ typeMismatches
  }

  private def validateRuntime(attributeExpression: WdlExpression, task: WdlTask, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Seq[SyntaxError] = {
    val invalidVariableReferences = for {
      variable <- referencesToAbsentValues(task, attributeExpression)
    } yield new SyntaxError(wdlSyntaxErrorFormatter.declarationContainsReferenceToAbsentValue(Option(task), variable))

    invalidVariableReferences.toSeq
  }

  private [wdl] def lookupType(from: Scope)(n: String): WomType = {
    val resolved = from.resolveVariable(n)
    resolved match {
      case Some(d: DeclarationInterface) => d.relativeWdlType(from)
      case Some(c: WdlCall) => WdlCallOutputsObjectType(c)
      case Some(s: Scatter) => s.collection.evaluateType(lookupType(s), new WdlStandardLibraryFunctionsType, Option(from)) match {
        case Success(WomArrayType(aType)) => aType
        // We don't need to check for a WOM map type, because
        // of the custom unapply in object WomArrayType
        case _ => throw new VariableLookupException(s"Variable $n references a scatter block ${s.fullyQualifiedName}, but the collection does not evaluate to an array")
      }
      case Some(_: WdlNamespace) => WdlNamespaceType
      case _ => throw new VariableLookupException(s"Could not resolve $n from scope ${from.fullyQualifiedName}")
    }
  }

  private def typeCheckDeclaration(decl: DeclarationInterface, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Option[SyntaxError] = {
    decl.expression flatMap { expr =>
      expr.evaluateType(lookupType(decl), new WdlStandardLibraryFunctionsType, Option(decl)) match {
        case Success(womType) =>
          if (!decl.womType.isCoerceableFrom(womType)) {
            Option(new SyntaxError(wdlSyntaxErrorFormatter.taskOutputExpressionTypeDoesNotMatchDeclaredType(
              declarationName(decl.ast), decl.womType, womType
            )))
          } else {
            expr.evaluate(NoLookup, NoFunctions) match {
              case Success(value) if decl.womType.coerceRawValue(value).isFailure =>
                Option(new SyntaxError(wdlSyntaxErrorFormatter.declarationExpressionNotCoerceableToTargetType(declarationName(decl.ast), decl.womType, value.womType)))
              case _ => None
            }
          }
        case Failure(_) => None
      }
    }
  }

  private def validateCallInputSection(call: WdlCall, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Seq[SyntaxError] = {
    val callInputSections = AstTools.callInputSectionIOMappings(call.ast, wdlSyntaxErrorFormatter)

    val invalidCallInputReferences = callInputSections flatMap { ast =>
      val lhs = ast.getAttribute("key").sourceString
      call.callable.inputNames.find(_ == lhs) match {
        case Some(_) => None
        case None => Option(new SyntaxError(wdlSyntaxErrorFormatter.callReferencesAbsentTaskInput(ast, call.callable.ast, lhs, call.unqualifiedName, call.isInstanceOf[WdlWorkflowCall])))
      }
    }

    def checkValidityOfMemberAccess(memberAccessAst: Ast): Option[SyntaxError] = {
      val memberAccess = MemberAccess(memberAccessAst)
      val requestedValue = memberAccess.rhsString

      // "Ignore local" so that we don't accidentally pick up an input assignment in the current call block that
      // happens to have the same name as the member access target (#3811). This probably nerfs ever fixing #4048.
      val resolvedScope: Option[Scope] = call.resolveVariable(memberAccess.lhsString, ignoreLocal = true)
      resolvedScope match {
        case Some(c: WdlCall) if c.outputs.exists(_.unqualifiedName == requestedValue) => None
        case Some(c: WdlCall) =>
          Option(new SyntaxError(wdlSyntaxErrorFormatter.memberAccessReferencesAbsentCallOutput(memberAccessAst, c)))
        case Some(s: Scatter) =>
          s.collection.evaluateType(lookupType(s), new WdlStandardLibraryFunctionsType) map {
            case WomArrayType(WomObjectType) => None
            case WomArrayType(_: WomPairType) if memberAccess.rhsString == "left" || memberAccess.rhsString == "right" => None
            // Maps get coerced into arrays of pairs, so this is also ok:
            case _: WomMapType if memberAccess.rhsString == "left" || memberAccess.rhsString == "right" => None
            case other => Option(new SyntaxError(wdlSyntaxErrorFormatter.badTargetTypeForMemberAccess(memberAccess, other)))
          } getOrElse None
        case Some(d: DeclarationInterface) => d.womType match {
          case _: WomPairType if memberAccess.rhsString == "left" || memberAccess.rhsString == "right" => None
          case WomObjectType => None
          case other => Option(new SyntaxError(wdlSyntaxErrorFormatter.badTargetTypeForMemberAccess(memberAccess, other)))
        }
        case Some(other) => Option(new SyntaxError(wdlSyntaxErrorFormatter.badTargetScopeForMemberAccess(memberAccess, other)))
        case None => None
          // In cases where there are many member accesses in the same Ast, it might be we can look up one layer but
          // not an inner layer.
          // It looks like we couldn't find this layer, so check whether there's an outer layer we can test access
          // for instead.
          memberAccess.lhsAst match {
            case outerMemberAccessAst: Ast if outerMemberAccessAst.getName == "MemberAccess" => checkValidityOfMemberAccess(outerMemberAccessAst)
            case _ => Option(new SyntaxError(wdlSyntaxErrorFormatter.noTargetForMemberAccess(memberAccess)))
          }
      }
    }

    /*
     * Ensures that the lhs corresponds to a call and the rhs corresponds to one of its outputs. We're only checking
     * top level MemberAccess ASTs because the sub-ASTs don't make sense w/o the context of the parent. For example
     * if we have "input: var=ns.ns1.my_task" it does not make sense to validate "ns1.my_task" by itself as it only
     * makes sense to validate that "ns.ns1.my_task" as a whole is coherent
     *
     */
    val invalidMemberAccesses: Seq[WdlParser.SyntaxError] = callInputSections flatMap { ast =>
      ast.getAttribute("value").findTopLevelMemberAccesses() flatMap checkValidityOfMemberAccess
    }

    invalidMemberAccesses ++
      invalidCallInputReferences
  }

  /**
    * Given a name, a collection of WdlNamespaces and a collection of Tasks, this method will attempt to find
    * a Task with that name within those collections.
    */
  def findTask(name: String, namespaces: Seq[WdlNamespace], tasks: Seq[WdlTask]): Option[WdlTask] = {
    findCallable(name, namespaces, tasks) collect { case t: WdlTask => t }
  }

  def findCallable(name: String, namespaces: Seq[WdlNamespace], callables: Seq[WdlCallable]): Option[WdlCallable] = {
    if (name.contains(".")) {
      val parts = name.split("\\.", 2)

      namespaces find (_.importedAs.contains(parts(0))) flatMap { x => findCallable(parts(1), x.namespaces, x.workflows ++ x.tasks) }
    } else callables.find(_.unqualifiedName == name)
  }


  private def readFile(wdlFile: Path): WorkflowSource = File(wdlFile).contentAsString

  def fileResolver(str: String): Draft2ResolvedImportBundle = {
    val path = Paths.get(str)
    Draft2ResolvedImportBundle(readFile(path), ResolvedImportRecord(path.toAbsolutePath.toString))
  }

  def directoryResolver(directory: File)(str: String): Draft2ResolvedImportBundle = {
    val absolutePathToFile = Paths.get(directory.path.resolve(str).toFile.getCanonicalPath)
    val absolutePathToImports = Paths.get(directory.toJava.getCanonicalPath)
    if (absolutePathToFile.startsWith(absolutePathToImports)) {
      fileResolver(absolutePathToFile.toString)
    } else {
      throw new IllegalArgumentException(s"$str is not a valid import")
    }
  }
}

object WdlNamespaceWithWorkflow {
  def load(workflowSource: WorkflowSource, importsResolvers: Seq[Draft2ImportResolver]): Try[WdlNamespaceWithWorkflow] = {
    from(WdlNamespace.loadUsingSource(workflowSource, None, Option(importsResolvers)))
  }

  @deprecated("To avoid unexpected default resolutions, I recommend using the load(String, Seq[ImportResolver] method of loading.", "23")
  def load(workflowSource: WorkflowSource): Try[WdlNamespaceWithWorkflow] = from(WdlNamespace.loadUsingSource(workflowSource, None, None))

  @deprecated("To avoid unexpected default resolutions, I recommend using the load(String, Seq[ImportResolver] method of loading.", "23")
  def load(workflowSource: WorkflowSource, importsDirectory: File): Try[WdlNamespaceWithWorkflow] = {
    val resolvers: Seq[Draft2ImportResolver] = Seq(WdlNamespace.directoryResolver(importsDirectory), WdlNamespace.fileResolver)
    load(workflowSource, resolvers)
  }

  @deprecated("To avoid unexpected default resolutions, I recommend using the load(String, Seq[ImportResolver] method of loading.", "23")
  def load(wdlFile: Path, importResolver: Draft2ImportResolver): Try[WdlNamespaceWithWorkflow] = from(WdlNamespace.loadUsingPath(wdlFile, None, Option(Seq(importResolver))))

  @deprecated("To avoid unexpected default resolutions, I recommend using the load(String, Seq[ImportResolver] method of loading.", "23")
  def load(workflowSource: WorkflowSource, importResolver: Draft2ImportResolver): Try[WdlNamespaceWithWorkflow] = {
    WdlNamespaceWithWorkflow.from(WdlNamespace.loadUsingSource(workflowSource, None, Option(Seq(importResolver))))
  }

  /**
    * Used to safely cast a WdlNamespace to a NamespaceWithWorkflow. Throws an IllegalArgumentException if another
    * form of WdlNamespace is passed to it
    */
  private def from(namespace: Try[WdlNamespace]): Try[WdlNamespaceWithWorkflow] = {
    namespace match {
      case Success(n: WdlNamespaceWithWorkflow) => Success(n)
      case Success(_) => Failure(new IllegalArgumentException("Namespace does not have a local workflow to run"))
      case Failure(f) => Failure(f)
    }
  }

  def apply(ast: Ast,
            workflow: WdlWorkflow,
            namespace: Option[String],
            imports: Seq[Import],
            namespaces: Seq[WdlNamespace],
            tasks: Seq[WdlTask],
            terminalMap: Map[Terminal, WorkflowSource],
            wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter,
            sourceString: String,
            importUri: Option[String],
            resolvedImportRecords: Set[ResolvedImportRecord]): WdlNamespaceWithWorkflow = {
    new WdlNamespaceWithWorkflow(namespace, workflow, imports, namespaces, tasks, terminalMap, wdlSyntaxErrorFormatter, ast, sourceString, importUri, resolvedImportRecords)
  }
}
