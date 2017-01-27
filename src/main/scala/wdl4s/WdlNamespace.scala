package wdl4s

import java.nio.file.{Path, Paths}

import better.files._
import lenthall.exception.AggregatedException
import lenthall.util.TryUtil
import wdl4s.AstTools.{AstNodeName, EnhancedAstNode}
import wdl4s.command.ParameterCommandPart
import wdl4s.exception._
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions, WdlStandardLibraryFunctionsType}
import wdl4s.parser.WdlParser._
import wdl4s.types._
import wdl4s.values._

import scala.collection.JavaConverters._
import scala.collection.mutable
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

/**
  * Represents a parsed WDL file
  */
sealed trait WdlNamespace extends WdlValue with Scope {
  final val wdlType = WdlNamespaceType
  def ast: Ast
  def resource = ast.findFirstTerminal.map(_.getResource).getOrElse("NONE")
  def importedAs: Option[String] // Used when imported with `as`
  def imports: Seq[Import]
  def namespaces: Seq[WdlNamespace]
  def tasks: Seq[Task]
  def workflows: Seq[Workflow]
  def terminalMap: Map[Terminal, WdlSource]
  def findTask(name: String): Option[Task] = tasks.find(_.name == name)
  override def unqualifiedName: LocallyQualifiedName = importedAs.getOrElse("")
  override def appearsInFqn: Boolean = importedAs.isDefined
  override def namespace: WdlNamespace = this
  def resolve(fqn: FullyQualifiedName): Option[Scope] = {
    (descendants + this).find(d => d.fullyQualifiedName == fqn || d.fullyQualifiedNameWithIndexScopes == fqn)
  }
  def resolveCallOrOutputOrDeclaration(fqn: FullyQualifiedName): Option[Scope] = {
    val callsAndOutputs = descendants collect { 
      case c: Call => c
      case d: Declaration => d
      case o: TaskOutput => o
      case o: CallOutput => o
    }
    callsAndOutputs.find(d => d.fullyQualifiedName == fqn || d.fullyQualifiedNameWithIndexScopes == fqn)
  }
}

/**
  * A WdlNamespace which doesn't have a locally defined Workflow.
  */
case class WdlNamespaceWithoutWorkflow(importedAs: Option[String],
                                       imports: Seq[Import],
                                       namespaces: Seq[WdlNamespace],
                                       tasks: Seq[Task],
                                       terminalMap: Map[Terminal, WdlSource],
                                       ast: Ast) extends WdlNamespace {
  val workflows = Seq.empty[Workflow]

}

/**
  * A WdlNamespace which has exactly one workflow defined.
  */
case class WdlNamespaceWithWorkflow(importedAs: Option[String],
                                    workflow: Workflow,
                                    imports: Seq[Import],
                                    namespaces: Seq[WdlNamespace],
                                    tasks: Seq[Task],
                                    terminalMap: Map[Terminal, WdlSource],
                                    wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter,
                                    ast: Ast) extends WdlNamespace {
  
  override val workflows = Seq(workflow)

  override def toString: String = s"[WdlNamespace importedAs=$importedAs]"


  /**
    * Confirm all required inputs are present and attempt to coerce raw inputs to `WdlValue`s.
    * This can fail if required raw inputs are missing or if the values for a specified raw input
    * cannot be coerced to the target type of the input as specified in the namespace.
    */
  def coerceRawInputs(rawInputs: WorkflowRawInputs): Try[WorkflowCoercedInputs] = {
    def coerceRawInput(input: WorkflowInput): Try[WdlValue] = input.fqn match {
      case _ if rawInputs.contains(input.fqn) =>
        val rawValue = rawInputs(input.fqn)
        input.wdlType.coerceRawValue(rawValue) match {
          case Success(value) => Success(value)
          case _ => Failure(new UnsatisfiedInputException(s"Could not coerce ${rawValue.getClass.getSimpleName} value for '${input.fqn}' ($rawValue) into: ${input.wdlType}"))
        }
      case _ =>
        input.optional match {
          case true => Success(WdlOptionalValue(input.wdlType.asInstanceOf[WdlOptionalType].memberType, None))
          case _ => Failure(new UnsatisfiedInputException(s"Required workflow input '${input.fqn}' not specified."))
        }
    }

    val tryCoercedValues = workflow.inputs map { case (fqn, input) => fqn -> coerceRawInput(input) }

    val (successes, failures) = tryCoercedValues.partition { case (_, tryValue) => tryValue.isSuccess }
    if (failures.isEmpty) {
      Try(for {
        (key, tryValue) <- successes
      } yield key -> tryValue.get)
    } else {
      val errors = failures.values.toList.collect { case f: Failure[_] => f.exception }
      Failure(ValidationException("Workflow input processing failed.", errors))
    }
  }

  /**
    * Some declarations need a value from the user and some have an expression attached to them.
    * For the declarations that have an expression attached to it already, evaluate the expression
    * and return the value. Only evaluates workflow level declarations. Other declarations will be evaluated at runtime.
    */
  def staticDeclarationsRecursive(userInputs: WorkflowCoercedInputs, wdlFunctions: WdlStandardLibraryFunctions): Try[WorkflowCoercedInputs] = {
    def evalDeclaration(accumulated: Map[FullyQualifiedName, Try[WdlValue]], current: Declaration): Map[FullyQualifiedName, Try[WdlValue]] = {
      current.expression match {
        case Some(expr) =>
          val successfulAccumulated = accumulated.collect({ case (k, v) if v.isSuccess => k -> v.get })
          val value = expr.evaluate(current.lookupFunction(successfulAccumulated ++ userInputs, wdlFunctions, NoOutputResolver, Map.empty[Scatter, Int]), wdlFunctions)
          val correctlyCoerced = value flatMap current.wdlType.coerceRawValue
          accumulated + (current.fullyQualifiedName -> correctlyCoerced)
        case None => accumulated
      }
    }

    def evalScope: Map[FullyQualifiedName, Try[WdlValue]] = {
      val workflowDeclarations = children.collect({ case w: Workflow => w.declarations }).flatten
      (declarations ++ workflowDeclarations).foldLeft(Map.empty[FullyQualifiedName, Try[WdlValue]])(evalDeclaration)
    }

    val filteredExceptions: Set[Class[_ <: Throwable]] = Set(classOf[OutputVariableLookupException], classOf[ScatterIndexNotFound])
    
    // Filter out declarations for which evaluation failed because a call output variable could not be resolved, or a shard could not be found,
    // as this method is meant for pre-execution validation
    val filtered = evalScope filterNot {
      case (_, Failure(ex)) if filteredExceptions.contains(ex.getClass) => true
      case (_, Failure(e: AggregatedException)) => e.throwables forall { ex => filteredExceptions.contains(ex.getClass) }
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
  * val namespace = WdlNamespace.load(new File("/path/to/file.wdl"))
  * namespace.workflow.calls foreach { call =>
  *   println(call)
  * }
  * }}}
  */
object WdlNamespace {

  def loadUsingPath(wdlFile: Path, resource: Option[String], importResolver: Option[Seq[ImportResolver]]): WdlNamespace = {
    load(readFile(wdlFile), resource.getOrElse(wdlFile.toString), importResolver.getOrElse(Seq(fileResolver)), None)
  }

  def loadUsingSource(wdlSource: WdlSource, resource: Option[String], importResolver: Option[Seq[ImportResolver]]): WdlNamespace = {
    load(wdlSource, resource.getOrElse("string"), importResolver.getOrElse(Seq(fileResolver)), None)
  }

  private def load(wdlSource: WdlSource, resource: String, importResolver: Seq[ImportResolver], importedAs: Option[String], root: Boolean = true): WdlNamespace = {
    WdlNamespace(AstTools.getAst(wdlSource, resource), wdlSource, importResolver, importedAs, root = root)
  }


  def apply(ast: Ast, source: WdlSource, importResolvers: Seq[ImportResolver], namespaceName: Option[String], root: Boolean = false): WdlNamespace = {

    val imports = for {
      importAst <- Option(ast).map(_.getAttribute("imports")).toSeq
      importStatement <- importAst.astListAsVector.map(Import(_))
    } yield importStatement

    def tryResolve(str: String, remainingResolvers: Seq[ImportResolver], errors: List[Throwable]): WdlSource = {
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

    /**
      * Translates all import statements to sub-namespaces
      */
    val namespaces: Seq[WdlNamespace] = for {
      imp <- imports
      source = tryResolve(imp.uri, importResolvers, List.empty)
    } yield WdlNamespace.load(source, imp.uri, importResolvers, Option(imp.namespaceName), root = false)

    /**
      * Map of Terminal -> WDL Source Code so the syntax error formatter can show line numbers
      */
    val terminalMap = AstTools.terminalMap(ast, source)
    val combinedTerminalMap = ((namespaces map { x => x.terminalMap }) ++ Seq(terminalMap)) reduce (_ ++ _)
    val wdlSyntaxErrorFormatter = WdlSyntaxErrorFormatter(combinedTerminalMap)

    val topLevelAsts = ast.getAttribute("body").astListAsVector.collect({ case a: Ast => a })

    /**
      * All `task` definitions of primary workflow.
      */
    val topLevelTasks: Seq[Task] = for {
      taskAst <- topLevelAsts if taskAst.getName == AstNodeName.Task
    } yield Task(taskAst, wdlSyntaxErrorFormatter)

    val workflows: Seq[Workflow] = for {
      workflowAst <- topLevelAsts if workflowAst.getName == AstNodeName.Workflow
    } yield Workflow(workflowAst, wdlSyntaxErrorFormatter)

    /**
      * Build scope tree recursively
      */
    val scopeIndexes: mutable.Map[Class[_ <: Scope], Int] = mutable.HashMap.empty.withDefaultValue(-1)

    def getScope(scopeAst: Ast, parent: Option[Scope]): Scope = {
      val scope = scopeAst.getName match {
        case AstNodeName.Call => Call(scopeAst, namespaces, topLevelTasks, workflows, wdlSyntaxErrorFormatter)
        case AstNodeName.Workflow => Workflow(scopeAst, wdlSyntaxErrorFormatter)
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
            case Some(task: Task) =>
              getScopeAsts(task.ast, "declarations").map(d => getScope(d, scope))
            case Some(workflow: Workflow) =>
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

    val namespace = workflows match {
      case Nil => WdlNamespaceWithoutWorkflow(namespaceName, imports, namespaces, topLevelTasks, terminalMap, ast)
      case Seq(w) => WdlNamespaceWithWorkflow(ast, w, namespaceName, imports, namespaces, topLevelTasks, terminalMap, wdlSyntaxErrorFormatter)

      case _ => throw new SyntaxError(wdlSyntaxErrorFormatter.tooManyWorkflows(ast.findAsts(AstNodeName.Workflow).asJava))
    }

    /**
      * Write-once var setting for parent/child relationships
      */

    def descendants(scope: Scope): Seq[Scope] = {
      val children = scope.children
      val childDescendants = scope.children.flatMap({
        case n: WdlNamespace => Seq.empty
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

    /**
      * SYNTAX CHECKS
      */

    val callInputSectionErrors = namespace.descendants.collect({ case c: TaskCall => c }).flatMap(
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

    val declarationErrors = for {
      descendant <- namespace.descendants
      declaration <- getDecls(descendant)
      error <- validateDeclaration(declaration, wdlSyntaxErrorFormatter)
    } yield error
    
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
      (namespace.descendants + namespace) collect { case workflow: Workflow => lookForDuplicates(workflow.outputs) }
    }

    val accumulatedErrors = expandedWorkflowOutputsDuplicationErrors ++ scopeDuplicationErrors

    val duplicateSiblingScopeNameErrors = accumulatedErrors.flatMap(_.errors).map(new SyntaxError(_)).toSeq

    val taskCommandReferenceErrors = for {
      task <- namespace.tasks
      param <- task.commandTemplate.collect({ case p: ParameterCommandPart => p })
      variable <- param.expression.variableReferences
      if !task.declarations.map(_.unqualifiedName).contains(variable.getSourceString)
    } yield new SyntaxError(wdlSyntaxErrorFormatter.commandExpressionContainsInvalidVariableReference(task.ast.getAttribute("name").asInstanceOf[Terminal], variable))

    val all = workflowOutputErrors ++ declarationErrors ++ callInputSectionErrors ++ taskCommandReferenceErrors ++ duplicateSiblingScopeNameErrors

    all.toSeq.sortWith({ case (l, r) => l.getMessage < r.getMessage }) match {
      case s: Seq[SyntaxError] if s.nonEmpty => throw s.head
      case _ =>
    }

    namespace
  }

  private def getDecls(scope: Scope): Seq[DeclarationInterface] = {
    scope match {
      case t: Task => t.declarations ++ t.outputs
      case w: Workflow => w.declarations ++ w.outputs
      case s => s.declarations
    }
  }

  private def declarationName(declarationAst: Ast): Terminal = declarationAst.getAttribute("name").asInstanceOf[Terminal]

  def validateDeclaration(declaration: DeclarationInterface, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Seq[SyntaxError] = {
    val invalidVariableReferences = for {
      expr <- declaration.expression.toSeq
      variable <- expr.variableReferences
      if declaration.resolveVariable(variable.sourceString).isEmpty
    } yield new SyntaxError(wdlSyntaxErrorFormatter.declarationContainsInvalidVariableReference(
      declarationName(declaration.ast),
      variable
    ))

    val typeMismatches = typeCheckDeclaration(declaration, wdlSyntaxErrorFormatter).toSeq

    invalidVariableReferences ++ typeMismatches
  }

  def lookupType(from: Scope)(n: String): WdlType = {
    val resolved = from.resolveVariable(n)
    resolved match {
      case Some(d: DeclarationInterface) => d.relativeWdlType(from)
      case Some(c: Call) => WdlCallOutputsObjectType(c)
      case Some(s: Scatter) => s.collection.evaluateType(lookupType(s), new WdlStandardLibraryFunctionsType, Option(from)) match {
        case Success(a: WdlArrayType) => a.memberType
        case _ => throw new VariableLookupException(s"Variable $n references a scatter block ${s.fullyQualifiedName}, but the collection does not evaluate to an array")
      }
      case Some(ns: WdlNamespace) => WdlNamespaceType
      case _ => throw new VariableLookupException(s"Could not resolve $n from scope ${from.fullyQualifiedName}")
    }
  }
  
  def typeCheckDeclaration(decl: DeclarationInterface, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Option[SyntaxError] = {
    decl.expression flatMap { expr =>
      expr.evaluateType(lookupType(decl), new WdlStandardLibraryFunctionsType, Option(decl)) match {
        case Success(wdlType) =>
          if (!decl.wdlType.isCoerceableFrom(wdlType)) {
            Option(new SyntaxError(wdlSyntaxErrorFormatter.taskOutputExpressionTypeDoesNotMatchDeclaredType(
              declarationName(decl.ast), decl.wdlType, wdlType
            )))
          } else {
            expr.evaluate(NoLookup, NoFunctions) match {
              case Success(value) if decl.wdlType.coerceRawValue(value).isFailure =>
                Option(new SyntaxError(wdlSyntaxErrorFormatter.declarationExpressionNotCoerceableToTargetType(
                  declarationName(decl.ast), decl.wdlType
                )))
              case _ => None
            }
          }
        case Failure(ex) => None
      }
    }
  }

  private def validateCallInputSection(call: TaskCall, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Seq[SyntaxError] = {
    val callInputSections = AstTools.callInputSectionIOMappings(call.ast, wdlSyntaxErrorFormatter)

    val invalidCallInputReferences = callInputSections flatMap { ast =>
      val lhs = ast.getAttribute("key").sourceString
      call.declarations.find(_.unqualifiedName == lhs) match {
        case Some(_) => None
        case None => Option(new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskInput(ast, call.task.ast)))
      }
    }

    /**
      * Ensures that the lhs corresponds to a call and the rhs corresponds to one of its outputs. We're only checking
      * top level MemberAccess ASTs because the sub-ASTs don't make sense w/o the context of the parent. For example
      * if we have "input: var=ns.ns1.my_task" it does not make sense to validate "ns1.my_task" by itself as it only
      * makes sense to validate that "ns.ns1.my_task" as a whole is coherent
      *
      * Run for its side effect (i.e. Exception) but we were previously using a Try and immediately calling .get on it
      * so it's the same thing
      */

    val invalidMemberAccesses = callInputSections flatMap { ast =>
      ast.getAttribute("value").findTopLevelMemberAccesses flatMap { memberAccessAst =>
        val memberAccess = MemberAccess(memberAccessAst)
        val resolvedScope: Option[Scope] = call.resolveVariable(memberAccess.lhs)
        resolvedScope match {
          case Some(c: Call) if c.outputs.exists(_.unqualifiedName == memberAccess.rhs) => None
          case Some(c: Call) =>
            Option(new SyntaxError(wdlSyntaxErrorFormatter.memberAccessReferencesBadCallInput(memberAccessAst, c)))
          case Some(s: Scatter) => 
            s.collection.evaluateType(lookupType(s), new WdlStandardLibraryFunctionsType) map {
              case WdlArrayType(WdlObjectType) => None
              case WdlArrayType(_: WdlPairType) if memberAccess.rhs == "left" || memberAccess.rhs == "right" => None
              case _ => Option(new SyntaxError(wdlSyntaxErrorFormatter.variableIsNotAnObject(memberAccessAst)))
            } getOrElse None
          case Some(d: Declaration) => d.wdlType match {
            case _: WdlPairType => None
            case _ => Option(new SyntaxError(wdlSyntaxErrorFormatter.variableIsNotAnObject(memberAccessAst)))
          }
          case None =>
            Option(new SyntaxError(wdlSyntaxErrorFormatter.undefinedMemberAccess(memberAccessAst)))
        }
      }
    }

    invalidMemberAccesses ++ invalidCallInputReferences
  }

  /**
    * Given a name, a collection of WdlNamespaces and a collection of Tasks, this method will attempt to find
    * a Task with that name within those collections.
    */
  def findTask(name: String, namespaces: Seq[WdlNamespace], tasks: Seq[Task]): Option[Task] = {
    findCallable(name, namespaces, tasks) collect { case t: Task => t }
  }

  def findCallable(name: String, namespaces: Seq[WdlNamespace], callables: Seq[Callable]): Option[Callable] = {
    if (name.contains(".")) {
      val parts = name.split("\\.", 2)

      namespaces find (_.importedAs.contains(parts(0))) flatMap { x => findCallable(parts(1), x.namespaces, x.workflows ++ x.tasks) }
    } else callables.find(_.unqualifiedName == name)
  }


  private def readFile(wdlFile: Path): WdlSource = File(wdlFile).contentAsString

  def fileResolver(str: String): WdlSource = readFile(Paths.get(str))

  def directoryResolver(directory: File)(str: String): WdlSource = {
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
  def load(wdlSource: WdlSource, importsResolvers: Seq[ImportResolver]): WdlNamespaceWithWorkflow = {
    from(WdlNamespace.loadUsingSource(wdlSource, None, Option(importsResolvers)))
  }

  @deprecated("To avoid unexpected default resolutions, I recommend using the load(String, Seq[ImportResolver] method of loading.", "23")
  def load(wdlSource: WdlSource): WdlNamespaceWithWorkflow = from(WdlNamespace.loadUsingSource(wdlSource, None, None))

  @deprecated("To avoid unexpected default resolutions, I recommend using the load(String, Seq[ImportResolver] method of loading.", "23")
  def load(wdlSource: WdlSource, importsDirectory: File): WdlNamespaceWithWorkflow = {
    val resolvers: Seq[ImportResolver] = Seq(WdlNamespace.directoryResolver(importsDirectory), WdlNamespace.fileResolver)
    load(wdlSource, resolvers)
  }

  @deprecated("To avoid unexpected default resolutions, I recommend using the load(String, Seq[ImportResolver] method of loading.", "23")
  def load(wdlFile: Path, importResolver: ImportResolver): WdlNamespaceWithWorkflow = from(WdlNamespace.loadUsingPath(wdlFile, None, Option(Seq(importResolver))))

  @deprecated("To avoid unexpected default resolutions, I recommend using the load(String, Seq[ImportResolver] method of loading.", "23")
  def load(wdlSource: WdlSource, importResolver: ImportResolver): WdlNamespaceWithWorkflow = {
    WdlNamespaceWithWorkflow.from(WdlNamespace.loadUsingSource(wdlSource, None, Option(Seq(importResolver))))
  }

  /**
    * Used to safely cast a WdlNamespace to a NamespaceWithWorkflow. Throws an IllegalArgumentException if another
    * form of WdlNamespace is passed to it
    */
  private def from(namespace: WdlNamespace): WdlNamespaceWithWorkflow = {
    namespace match {
      case n: WdlNamespaceWithWorkflow => n
      case _ => throw new IllegalArgumentException("Namespace does not have a local workflow to run")
    }
  }

  def apply(ast: Ast, workflow: Workflow, namespace: Option[String], imports: Seq[Import],
            namespaces: Seq[WdlNamespace], tasks: Seq[Task], terminalMap: Map[Terminal, WdlSource],
            wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): WdlNamespaceWithWorkflow = {
    new WdlNamespaceWithWorkflow(namespace, workflow, imports, namespaces, tasks, terminalMap, wdlSyntaxErrorFormatter, ast)
  }
}
