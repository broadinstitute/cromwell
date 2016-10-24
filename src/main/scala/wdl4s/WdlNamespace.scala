package wdl4s

import java.nio.file.{Path, Paths}

import better.files._
import cats.data.NonEmptyList
import wdl4s.AstTools.{AstNodeName, EnhancedAstNode}
import wdl4s.command.ParameterCommandPart
import wdl4s.expression.{NoFunctions, WdlStandardLibraryFunctions, WdlStandardLibraryFunctionsType}
import wdl4s.parser.WdlParser._
import wdl4s.types._
import wdl4s.util.TryUtil
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
  val workflows = Seq(workflow)

  override def toString: String = s"[WdlNamespace importedAs=$importedAs]"

  /**
    * Confirm all required inputs are present and attempt to coerce raw inputs to `WdlValue`s.
    * This can fail if required raw inputs are missing or if the values for a specified raw input
    * cannot be coerced to the target type of the input as specified in the namespace.
    */
  def coerceRawInputs(rawInputs: WorkflowRawInputs): Try[WorkflowCoercedInputs] = {
    def coerceRawInput(input: WorkflowInput): Try[Option[WdlValue]] = input.fqn match {
      case _ if rawInputs.contains(input.fqn) =>
        val rawValue = rawInputs(input.fqn)
        input.wdlType.coerceRawValue(rawValue) match {
          case Success(value) => Success(Some(value))
          case _ => Failure(new UnsatisfiedInputsException(s"Could not coerce value for '${input.fqn}' into: ${input.wdlType}"))
        }
      case _ =>
        input.optional match {
          case true => Success(None)
          case _ => Failure(new UnsatisfiedInputsException(s"Required workflow input '${input.fqn}' not specified."))
        }
    }

    val tryCoercedValues = workflow.inputs map { case (fqn, input) => fqn -> coerceRawInput(input) }

    val (successes, failures) = tryCoercedValues.partition { case (_, tryValue) => tryValue.isSuccess }
    if (failures.isEmpty) {
      Try(for {
        (key, tryValue) <- successes
        optionValue = tryValue.get if tryValue.get.isDefined
      } yield key -> optionValue.get)
    } else {
      val errors = failures.values.collect { case f: Failure[_] => f.exception.getMessage }
      // .fromListUnsafe because failures is guaranteed to be nonEmpty
      Failure(new ValidationException("Workflow input processing failed.", NonEmptyList.fromListUnsafe(errors.toList)))
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
          accumulated + (current.fullyQualifiedName -> value)
        case None => accumulated
      }
    }

    def evalScope(scope: Scope): Map[FullyQualifiedName, Try[WdlValue]] = {
      val evaledDeclarations = scope.declarations.foldLeft(Map.empty[FullyQualifiedName, Try[WdlValue]])(evalDeclaration)
      val nonTasks = scope.children.collect({
        case w: Workflow => w
      })
      evaledDeclarations ++ nonTasks.flatMap(evalScope).toMap
    }

    TryUtil.sequenceMap(evalScope(this))
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
  def load(wdlFile: Path): WdlNamespace = {
    load(readFile(wdlFile), wdlFile.toString, localImportResolver, None)
  }

  def load(wdlFile: Path, importResolver: ImportResolver): WdlNamespace = {
    load(readFile(wdlFile), wdlFile.toString, importResolver, None)
  }

  def load(wdlSource: WdlSource): WdlNamespace = {
    load(wdlSource, "string", localImportResolver, None)
  }

  def load(wdlSource: WdlSource, importResolver: ImportResolver): WdlNamespace = {
    load(wdlSource, "string", importResolver, None)
  }

  def load(wdlSource: WdlSource, resource: String): WdlNamespace = {
    load(wdlSource, resource, localImportResolver, None)
  }

  def load(wdlSource: WdlSource, resource: String, importResolver: ImportResolver): WdlNamespace = {
    load(wdlSource, resource, importResolver, None)
  }

  private def load(wdlSource: WdlSource, resource: String, importResolver: ImportResolver, importedAs: Option[String]): WdlNamespace = {
    WdlNamespace(AstTools.getAst(wdlSource, resource), wdlSource, importResolver, importedAs, root = true)
  }

  def apply(ast: Ast, source: WdlSource, importResolver: ImportResolver, namespaceName: Option[String], root: Boolean = false): WdlNamespace = {
    val imports = for {
      importAst <- Option(ast).map(_.getAttribute("imports")).toSeq
      importStatement <- importAst.astListAsVector.map(Import(_))
    } yield importStatement

    /**
      * For import statements with an 'as' clause, this translates to a sub-namespace
      */
    val namespaces: Seq[WdlNamespace] = for {
      imp <- imports if imp.namespaceName.isDefined
      source = importResolver(imp.uri)
      if source.length > 0
    } yield WdlNamespace.load(source, imp.uri, importResolver, imp.namespaceName)

    /**
      * For import statements without an 'as' clause, this is treated as a c-style #include
      */
    val nonNamespacedImports = for {
      imp <- imports if imp.namespaceName.isEmpty
      source = importResolver(imp.uri)
      ast = AstTools.getAst(source, imp.uri)
      scopeAst <- ast.getAttribute("body").astListAsVector.collect({ case a: Ast => a })
    } yield ast -> source

    val topLevelScopeAsts = for {
      namespaceAst <- nonNamespacedImports.map(_._1) :+ ast
      scopeAst <- namespaceAst.getAttribute("body").astListAsVector.collect({ case a: Ast => a })
    } yield scopeAst

    /**
      * Map of Terminal -> WDL Source Code so the syntax error formatter can show line numbers
      */
    val nonNamespacedImportsTerminalMap = nonNamespacedImports.flatMap({ case (a, s) => AstTools.terminalMap(a, s) }).toMap
    val terminalMap = nonNamespacedImportsTerminalMap ++ AstTools.terminalMap(ast, source)
    val combinedTerminalMap = ((namespaces map {x => x.terminalMap}) ++ Seq(terminalMap)) reduce (_ ++ _)
    val wdlSyntaxErrorFormatter = new WdlSyntaxErrorFormatter(combinedTerminalMap)

    /**
      * All `task` definitions, including ones from import statements with no 'as' clause
      */
    val tasks: Seq[Task] = for {
      taskAst <- topLevelScopeAsts if taskAst.getName == AstNodeName.Task
    } yield Task(taskAst, wdlSyntaxErrorFormatter)

    /**
      * Build scope tree recursively
      */
    val scopeIndexes: mutable.Map[Class[_ <: Scope], Int] = mutable.HashMap.empty.withDefaultValue(-1)

    def getScope(scopeAst: Ast, parent: Option[Scope]): Scope = {
      val scope = scopeAst.getName match {
        case AstNodeName.Call => Call(scopeAst, namespaces, tasks, wdlSyntaxErrorFormatter)
        case AstNodeName.Workflow => Workflow(scopeAst, wdlSyntaxErrorFormatter)
        case AstNodeName.Declaration => Declaration(scopeAst, wdlSyntaxErrorFormatter, parent)
        case AstNodeName.Scatter =>
          scopeIndexes(classOf[Scatter]) += 1
          Scatter(scopeAst, scopeIndexes(classOf[Scatter]))
        case AstNodeName.If =>
          scopeIndexes(classOf[If]) += 1
          If(scopeAst, scopeIndexes(classOf[If]))
        case AstNodeName.Output => TaskOutput(scopeAst, wdlSyntaxErrorFormatter, parent)
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

      scopeAst.getName match {
        case AstNodeName.Task => getTaskInputsOutputs(scopeAst)
        case AstNodeName.Declaration | AstNodeName.Output => Seq.empty[Scope]
        case AstNodeName.Call =>
          val referencedTask = findTask(scopeAst.getAttribute("task").sourceString, namespaces, tasks)
          referencedTask match {
            case Some(task) => getScopeAsts(task.ast, "declarations").map(d => getScope(d, scope))
            case None => Seq.empty[Scope]
          }
        case AstNodeName.Workflow | AstNodeName.Scatter | AstNodeName.If | AstNodeName.Namespace =>
          getScopeAsts(scopeAst, "body").map(getScope(_, scope))
      }
    }

    val nonTaskScopes = for {
      ast <- topLevelScopeAsts
      if ast.getName != AstNodeName.Task
    } yield ast

    val children = tasks ++ namespaces ++ nonTaskScopes.map(ast => getScope(ast, parent = None))

    val namespace = children.collect({ case w: Workflow => w }) match {
      case Nil => WdlNamespaceWithoutWorkflow(namespaceName, imports, namespaces, tasks, terminalMap, ast)
      case Seq(workflow) => WdlNamespaceWithWorkflow(ast, workflow, namespaceName, imports, namespaces, tasks, terminalMap, wdlSyntaxErrorFormatter)
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

    tasks foreach { task =>
      task.children = getChildren(task.ast, Option(task))
      task.children.foreach(_.parent = task)
    }

    descendants(namespace).foreach(_.namespace = namespace)

    /**
      * SYNTAX CHECKS
      */

    val callInputSectionErrors = namespace.descendants.collect({ case c: Call => c }).flatMap(
      validateCallInputSection(_, wdlSyntaxErrorFormatter)
    )

    val declarationErrors = for {
      descendant <- namespace.descendants
      declaration <- getDecls(descendant)
      error <- validateDeclaration(declaration, wdlSyntaxErrorFormatter)
    } yield error

    def scopeNameAndTerminal(scope: Scope): (String, Terminal) = {
      scope match {
        case ns: WdlNamespace => ("Namespace", imports.find(_.namespaceName == ns.importedAs).get.namespaceTerminal.get)
        case s: Scope => (s.getClass.getSimpleName, s.ast.findFirstTerminal.get)
      }
    }
    case class ScopeAccumulator(accumulated: Seq[Scope] = Seq.empty, errors: Seq[String] = Seq.empty)
    val accumulatedErrors = (namespace.descendants + namespace).filter(_.namespace == namespace) map { scope =>

      scope.children.foldLeft(ScopeAccumulator()) { (acc, cur) =>
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
    val duplicateSiblingScopeNameErrors = accumulatedErrors.flatMap(_.errors).map(new SyntaxError(_)).toSeq

    val taskCommandReferenceErrors = for {
      task <- namespace.tasks
      param <- task.commandTemplate.collect({ case p: ParameterCommandPart => p })
      variable <- param.expression.variableReferences
      if !task.declarations.map(_.unqualifiedName).contains(variable.getSourceString)
    } yield new SyntaxError(wdlSyntaxErrorFormatter.commandExpressionContainsInvalidVariableReference(task.ast.getAttribute("name").asInstanceOf[Terminal], variable))

    val all = declarationErrors ++ callInputSectionErrors ++ taskCommandReferenceErrors ++ duplicateSiblingScopeNameErrors

    all.toSeq.sortWith({ case (l, r) => l.getMessage < r.getMessage}) match {
      case s: Seq[SyntaxError] if s.nonEmpty => throw s.head
      case _ =>
    }

    namespace
  }

  private def getDecls(scope: Scope): Seq[DeclarationInterface] = {
    scope match {
      case t: Task => t.declarations ++ t.outputs
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

  def typeCheckDeclaration(decl: DeclarationInterface, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Option[SyntaxError] = {
    def lookupType(from: Scope)(n: String): WdlType = {
      from.resolveVariable(n) match {
        case Some(d: DeclarationInterface) => d.wdlType
        case Some(c: Call) => WdlObjectType
        case Some(s: Scatter) => s.collection.evaluateType(lookupType(s), new WdlStandardLibraryFunctionsType) match {
          case Success(a: WdlArrayType) => a.memberType
          case _ => throw new VariableLookupException(s"Variable $n references a scatter block ${s.fullyQualifiedName}, but the collection does not evaluate to an array")
        }
        case Some(ns: WdlNamespace) => WdlNamespaceType
        case _ => throw new VariableLookupException(s"Could not resolve $n from scope ${from.fullyQualifiedName}")
      }
    }

    decl.expression flatMap { expr =>
      expr.evaluateType(lookupType(decl), new WdlStandardLibraryFunctionsType) match {
        case Success(wdlType) if !decl.wdlType.isCoerceableFrom(wdlType) =>
          Option(new SyntaxError(wdlSyntaxErrorFormatter.taskOutputExpressionTypeDoesNotMatchDeclaredType(
            declarationName(decl.ast), wdlType, decl.wdlType
          )))
        case Success(wdlType) =>
          expr.evaluate(NoLookup, NoFunctions) match {
            case Success(value) if decl.wdlType.coerceRawValue(value).isFailure =>
              Option(new SyntaxError(wdlSyntaxErrorFormatter.declarationExpressionNotCoerceableToTargetType(
                declarationName(decl.ast), decl.wdlType
              )))
            case _ => None
          }
        case Failure(ex) => None
      }
    }
  }

  private def validateCallInputSection(call: Call, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Seq[SyntaxError] = {
    val callInputSections = AstTools.callInputSectionIOMappings(call.ast, wdlSyntaxErrorFormatter)

    val invalidCallInputReferences = callInputSections flatMap { ast =>
      val lhs = ast.getAttribute("key").sourceString
      val rhs = ast.getAttribute("value")
      call.declarations.find(_.unqualifiedName == lhs) match {
        case Some(decl) => None
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
        call.resolveVariable(memberAccess.lhs) match {
          case Some(c: Call) if c.task.outputs.exists(_.unqualifiedName == memberAccess.rhs) => None
          case Some(c: Call) =>
            Option(new SyntaxError(wdlSyntaxErrorFormatter.memberAccessReferencesBadTaskInput(memberAccessAst)))
          case None =>
            Option(new SyntaxError(wdlSyntaxErrorFormatter.undefinedMemberAccess(memberAccessAst)))
        }
      }
    }

    invalidMemberAccesses ++ invalidCallInputReferences
  }

  /**
    * Given a name, a collection of WdlNamespaces and a collection of Tasks will attempt to find a Task
    * with that name within the WdlNamespaces
    */
  def findTask(name: String, namespaces: Seq[WdlNamespace], tasks: Seq[Task]): Option[Task] = {
    if (name.contains(".")) {
      val parts = name.split("\\.", 2)
      /* This is supposed to resolve a dot-notation'd string (e.g. "a.b.c") by recursively
       * traversing child namespaces or resolving to a task.
       *
       * For example:
       * findTasks("a.b.c") would first find namespace "a" and then return a.findTasks("b.c")
       * a.findTasks("b.c") would call a.b.findTasks("c")
       * a.b.findTasks("c") would return the task named "c" in the "b" namespace
       */
      namespaces find (_.importedAs.contains(parts(0))) flatMap { x => findTask(parts(1), x.namespaces, x.tasks)}
    } else tasks.find(_.name == name)
  }

  private def localImportResolver(path: String): WdlSource = readFile(Paths.get(path))

  private def readFile(wdlFile: Path): WdlSource = File(wdlFile).contentAsString
}

object WdlNamespaceWithWorkflow {
  def load(wdlSource: WdlSource): WdlNamespaceWithWorkflow = from(WdlNamespace.load(wdlSource))

  def load(wdlFile: Path, importResolver: ImportResolver): WdlNamespaceWithWorkflow = from(WdlNamespace.load(wdlFile, importResolver))

  def load(wdlSource: WdlSource, importResolver: ImportResolver): WdlNamespaceWithWorkflow = {
    WdlNamespaceWithWorkflow.from(WdlNamespace.load(wdlSource, importResolver))
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
