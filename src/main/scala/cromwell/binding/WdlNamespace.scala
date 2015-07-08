package cromwell.binding

import java.io.File

import cromwell.binding.WdlNamespace.ImportResolver
import cromwell.binding.command._
import cromwell.binding.types._
import cromwell.binding.values.WdlValue
import cromwell.parser.WdlParser
import cromwell.parser.WdlParser._
import cromwell.util.FileUtil
import cromwell.binding.AstTools.EnhancedAstNode
import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}
import scala.language.postfixOps

/**
 * Main interface into the `cromwell.binding` package.
 *
 * Example usage:
 *
 * {{{
 * val namespace = WdlNamespace.process(new File("/path/to/file.wdl"))
 * binding.workflow.calls foreach { call =>
 *      println(call)
 * }
 * }}}
 */

object WdlNamespace {
  type ImportResolver = String => WdlSource
  def localImportResolver(path: String): WdlSource = readFile(new File(path))

  /**
   * Given a pointer to a WDL file, parse the text and build Workflow and Task
   * objects.
   *
   * @param wdlFile The file to parse/process
   * @return WdlBinding object with the parsed results
   * @throws WdlParser.SyntaxError if there was a problem parsing the source code
   * @throws UnsupportedOperationException if an error occurred constructing the
   *                                       Workflow and Task objects
   *
   */
  def load(wdlFile: File): WdlNamespace =
    load(readFile(wdlFile), wdlFile.toString, localImportResolver, None)

  def load(wdlFile: File, importResolver: ImportResolver): WdlNamespace =
    load(readFile(wdlFile), wdlFile.toString, importResolver, None)

  def load(wdlSource: WdlSource): WdlNamespace =
    load(wdlSource, "string", localImportResolver, None)

  def load(wdlSource: WdlSource, importResolver: ImportResolver): WdlNamespace =
    load(wdlSource, "string", importResolver, None)

  def load(wdlSource: WdlSource, resource: String): WdlNamespace =
    load(wdlSource, resource, localImportResolver, None)

  def load(wdlSource: WdlSource, resource: String, importResolver: ImportResolver): WdlNamespace =
    load(wdlSource, resource, importResolver, None)

  private def load(wdlSource: WdlSource, resource: String, importResolver: ImportResolver, namespace: Option[String]): WdlNamespace =
    new WdlNamespace(AstTools.getAst(wdlSource, resource), wdlSource, importResolver, namespace)

  def readFile(wdlFile: File): WdlSource = FileUtil.slurp(wdlFile)
}

case class WdlNamespace(ast: Ast, source: WdlSource, importResolver: ImportResolver, namespace: Option[String]) extends WdlValue {
  val wdlType = WdlNamespaceType

  import AstTools.AstNodeName

  val invalidTask = Task("INVALID", Command(Seq.empty[CommandPart]), Seq.empty[TaskOutput], Map.empty[String, String])

  /**
   * All `import` statement strings at the top of the document
   */
  val imports = ast.getAttribute("imports").asInstanceOf[AstList].asScala.toVector.map { i =>
    val uri = i.asInstanceOf[Ast].getAttribute("uri").sourceString()
    val namespaceAst = i.asInstanceOf[Ast].getAttribute("namespace")
    val namespace = if (namespaceAst == null) None else Option(namespaceAst.sourceString())
    Import(uri, namespace)
  }

  /* WdlBinding objects for each import statement */
  val namespaces = for {
    i <- imports
    source = importResolver(i.uri)
    if source.length > 0
  } yield WdlNamespace.load(source, i.uri, importResolver, i.namespace)

  /* Create a map of Terminal -> WdlBinding */
  val terminalMap = AstTools.terminalMap(ast, source)
  val combinedTerminalMap = ((namespaces map {x => x.terminalMap}) ++ Seq(terminalMap)) reduce (_ ++ _)
  val wdlSyntaxErrorFormatter = new WdlSyntaxErrorFormatter(combinedTerminalMap)

  /**
   * All imported `task` definitions for `import` statements without a namespace (e.g. no `as` clause)
   * These tasks are considered to be in this current workspace
   */
  val importedTaskAsts: Seq[Ast] = namespaces flatMap { b =>
    b.namespace match {
      case None => b.taskAsts
      case _ => Seq.empty[Ast]
    }
  }
  val importedTasks: Seq[Task] = namespaces flatMap { b =>
    b.namespace match {
      case None => b.tasks
      case _ => Seq.empty[Task]
    }
  }

  /**
   * All `task` definitions defined in the WDL file (i.e. not imported)
   */
  val localTaskAsts = ast.findAsts(AstNodeName.Task)
  val localTasks = localTaskAsts.map(processTask)

  /**
   * All `task` definitions, including local and imported ones
   */
  val taskAsts = localTaskAsts ++ importedTaskAsts
  val tasks = localTasks ++ importedTasks

  /**
   * All imported `Workflow`s
   */
  val importedWorkflows: Seq[Workflow] = namespaces flatMap { b =>
    b.namespace match {
      case None => b.workflows
      case _ => Seq.empty[Workflow]
    }
  }

  /**
   * All `Workflow`s defined in the WDL file (i.e. not imported)
   */
  val localWorkflows = ast.findAsts(AstNodeName.Workflow).map(processWorkflow)

  /**
   * All `Workflow` definitions, including local and imported ones
   */
  val workflows = importedWorkflows ++ localWorkflows

  lazy val workflow = workflows.head

  lazy val calls = workflow.calls

  /* All `Task`s and `Workflow`s */
  val executables = workflows ++ tasks

  workflows foreach { workflow =>
    workflow.calls.foreach {
      _.setParent(workflow)
    }
  }

  def findTaskAst(name: String): Option[Ast] = {
    if (name.contains(".")) {
      val parts = name.split("\\.", 2)
      namespaces find {_.namespace == Some(parts(0))} flatMap {_.findTaskAst(parts(1))}
    } else {
      taskAsts.find{t => t.getAttribute("name").sourceString == name}
    }
  }

  def findTask(name: String): Option[Task] = {
    if (name.contains(".")) {
      val parts = name.split("\\.", 2)
      namespaces find {_.namespace == Some(parts(0))} flatMap {_.findTask(parts(1))}
    } else {
      tasks.find(_.name == name)
    }
  }

  private def processWorkflow(ast: Ast): Workflow = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val calls = ast.findAsts(AstNodeName.Call) map processCall
    new Workflow(name, calls)
  }

  private def processTask(ast: Ast): Task = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val commandAsts = ast.findAsts(AstNodeName.Command)
    if (commandAsts.size != 1) throw new UnsupportedOperationException("Expecting only one Command AST")
    val command = processCommand(commandAsts.head.getAttribute("parts").asInstanceOf[AstList])
    val outputs = ast.findAsts(AstNodeName.Output) map processTaskOutput
    new Task(name, command, outputs, buildRuntimeAttributes(ast))
  }

  private def buildRuntimeAttributes(ast: Ast): RuntimeAttributes = {
    val asts = ast.findAsts(AstNodeName.Runtime)
    if (asts.size > 1) throw new UnsupportedOperationException("Only one runtime block may be defined per task")
    val astList = asts.headOption map {_.getAttribute("map").asInstanceOf[AstList]}
    astList map processRuntimeAttributes getOrElse Map.empty[String, String]
  }

  private def processRuntimeAttributes(astList: AstList): RuntimeAttributes = {
    astList.asScala.toVector map {a => processRuntimeAttribute(a.asInstanceOf[Ast])} toMap
  }

  private def processRuntimeAttribute(ast: Ast): RuntimeAttribute = {
    (ast.getAttribute("key").sourceString(), ast.getAttribute("value").sourceString())
  }

  private def processCommand(astList: AstList): Command = {
    val parts = astList.asScala.toVector.map {
      case x: Terminal => new StringCommandPart(x.getSourceString)
      case x: Ast => processCommandParameter(x)
    }
    new Command(parts)
  }

  private def processCommandParameter(ast: Ast): ParameterCommandPart = {
    val wdlType = processWdlType(ast.getAttribute("type"))
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    new ParameterCommandPart(wdlType, name)
  }

  private def processTaskOutput(ast: Ast): TaskOutput = {
    val wdlType = processWdlType(ast.getAttribute("type"))
    val name = ast.getAttribute("var").sourceString()
    val expression = ast.getAttribute("expression")
    new TaskOutput(name, wdlType, new WdlExpression(expression))
  }

  private def processCall(ast: Ast): Call = {
    val alias: Option[String] = ast.getAttribute("alias") match {
      case x: Terminal => Option(x.getSourceString)
      case _ => None
    }
    val taskName = ast.getAttribute("task").sourceString()
    val task = findTask(taskName) getOrElse invalidTask
    val inputs = processCallInput(ast)
    new Call(alias, taskName, task, inputs.toMap, this)
  }

  private def processCallInput(ast: Ast): Map[String, WdlExpression] = AstTools.callInputAsts(ast).map {a =>
    val key = a.getAttribute("key").sourceString()
    val expression = new WdlExpression(a.getAttribute("value"))
    (key, expression)
  }.toMap

  private def processWdlType(ast: AstNode): WdlType = {
    ast match {
      case t: Terminal =>
        t.getSourceString match {
          case WdlFileType.toWdlString => WdlFileType
          case WdlStringType.toWdlString => WdlStringType
          case WdlIntegerType.toWdlString => WdlIntegerType
          case WdlFloatType.toWdlString => WdlFloatType
          case WdlBooleanType.toWdlString => WdlBooleanType
          case WdlObjectType.toWdlString => WdlObjectType
        }
      case null => WdlStringType
      case _ => throw new UnsupportedOperationException("Implement this later for compound types")
    }
  }

  /**
   * Confirm all required inputs are present and attempt to coerce raw inputs to `WdlValue`s.
   * This can fail if required raw inputs are missing or if the values for a specified raw input
   * cannot be coerced to the target type of the input as specified in the namespace.
   */
  def coerceRawInputs(rawInputs: WorkflowRawInputs): Try[WorkflowCoercedInputs] = {

    def coerceRawInput(fqn: FullyQualifiedName, wdlType: WdlType): Try[WdlValue] = fqn match {
      case _ if rawInputs.contains(fqn) =>
        val rawInput = rawInputs.get(fqn).get
        wdlType.coerceRawValue(rawInput).recoverWith {
          case e => Failure(new UnsatisfiedInputsException(s"Failed to coerce input $fqn value $rawInput of ${rawInput.getClass} to $wdlType."))
        }
      case _ => Failure(new UnsatisfiedInputsException(s"Required workflow input '$fqn' not specified."))
    }

    val tryCoercedValues = workflows.head.inputs.map { case (fqn, wdlType) =>
      fqn -> coerceRawInput(fqn, wdlType)
    }

    val (successes, failures) = tryCoercedValues.partition { case (_, tryValue) => tryValue.isSuccess }
    if (failures.isEmpty) {
      Try(successes.map { case (key, tryValue) => key -> tryValue.get })
    } else {
      val message = failures.values.collect { case f: Failure[_] => f.exception.getMessage }.mkString("\n")
      Failure(new UnsatisfiedInputsException(message))
    }
  }

  /**
   * Validate the following things about the AST:
   *
   * 1) All `Call` blocks reference tasks that exist
   * 2) All `Call` inputs reference actual variables on the corresponding task
   * 3) Tasks do not have duplicate inputs
   * 4) Calls do not reference the same task input more than once
   * 5) Tasks in this namespace have unique names
   * 6) Tasks and namespaces don't have overlapping names
   * 7) `Call` input expressions (right-hand side) should only use the MemberAccess
   * syntax (e.g: x.y) on WdlObjects (which include other `Call` invocations)
   * 8) `Call` input expressions (right-hand side) should only reference identifiers
   * that will resolve when evaluated
   */
  private def validate(): Unit = {
    val workflowAsts = ast.findAsts(AstNodeName.Workflow)

    if (workflowAsts.size > 1) {
      throw new SyntaxError(wdlSyntaxErrorFormatter.tooManyWorkflows(workflowAsts.asJava))
    }

    tasks foreach {task =>
      val taskAstsWithSameName = taskAsts filter {_.getAttribute("name").sourceString == task.name}
      /* No two tasks can have the same name */
      if (taskAstsWithSameName.size > 1) {
        throw new SyntaxError(wdlSyntaxErrorFormatter.duplicateTask(taskAstsWithSameName))
      }
      /* A task can not have duplicated inputs */
      val commandLineInputs = taskAstsWithSameName.head.findAsts(AstNodeName.CommandParameter)
      commandLineInputs foreach {input =>
        val inputName = input.getAttribute("name").sourceString()
        val inputsWithSameName = commandLineInputs filter {_.getAttribute("name").sourceString == inputName}
        if (inputsWithSameName.size > 1) {
          throw new SyntaxError(wdlSyntaxErrorFormatter.taskHasDuplicatedInputs(taskAstsWithSameName.head, inputsWithSameName))
        }
      }
    }

    /* Ensure that no namespaces collide with task/workflow names */
    ast.getAttribute("imports").asInstanceOf[AstList].asScala.toVector.foreach {i =>
      val namespaceAst = i.asInstanceOf[Ast].getAttribute("namespace").asInstanceOf[Terminal]
      if (namespaceAst != null) {
        findTaskAst(namespaceAst.sourceString()) match {
          case Some(taskAst) =>
            throw new SyntaxError(wdlSyntaxErrorFormatter.taskAndNamespaceHaveSameName(taskAst, namespaceAst))
          case _ =>
        }
        workflowAsts foreach {workflowAst =>
          if (namespaceAst.sourceString == workflowAst.getAttribute("name").sourceString) {
            throw new SyntaxError(wdlSyntaxErrorFormatter.workflowAndNamespaceHaveSameName(workflowAst, namespaceAst))
          }
        }
      }
    }

    workflowAsts foreach { workflowAst =>
      val workflow = localWorkflows.head
      workflowAst.findAsts(AstNodeName.Call) foreach { callAst =>
        val taskName = callAst.getAttribute("task").sourceString()
        val taskAst = findTaskAst(taskName) getOrElse {
          throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskName(callAst, taskName))
        }
        val task = findTask(taskName) getOrElse {
          throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskName(callAst, taskName))
        }

        val callName = Option(callAst.getAttribute("alias")).getOrElse(callAst.getAttribute("task")).sourceString()

        val call = workflow.calls.find{_.name == callName} getOrElse {
          throw new SyntaxError(s"Cannot find call with name $callName")
        }

        call.inputMappings foreach { inputKv =>
          task.inputs find { taskInput => taskInput._1 == inputKv._1 } getOrElse {
            val callInput = AstTools.callInputAsts(callAst) find {
              _.getAttribute("key").sourceString == inputKv._1
            } getOrElse {
              throw new SyntaxError(s"Can't find call input: ${inputKv._1}")
            }
            throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskInput(callInput, taskAst))
          }

          /* All MemberAccess ASTs that are not contained in other MemberAccess ASTs */
          inputKv._2.ast.findTopLevelMemberAccesses().map(getCallFromMemberAccessAst).map {_.get}
        }
      }
    }
  }

  /* Partially evaluate MemberAccess ASTs to make sure they're not nonsense at compile time */
  def getCallFromMemberAccessAst(ast: Ast): Try[Call] = {
    def callFromName(name: String): Try[Call] = {
      workflows.head.calls.find(_.name == name) match {
        case Some(c:Call) => Success(c)
        case _ => Failure(new SyntaxError(wdlSyntaxErrorFormatter.undefinedMemberAccess(ast)))
      }
    }
    val rhs = ast.getAttribute("rhs").sourceString()

    /**
     * The left-hand side of a member-access AST should always be interpreted as a String
     * Sometimes, the left-hand side is itself a MemberAccess AST, like in the expression
     * for `call t1` below.  In that example, callFromName("ns.ns2.task_name") would be
     * called.  In the `call t2` example, callFromName("alias") is called
     *
     * import "test.wdl" as ns
     * workflow w {
     *  call ns.ns2.task_name
     *  call t1 {
     *    input: x=ns.ns2.task_name.output
     *  }
     *
     *  call ns.ns2.task_name as alias
     *  call t2 {
     *    input: y=alias.output
     *  }
     *}
     */
    val lhs = callFromName(ast.getAttribute("lhs") match {
      case a: Ast => WdlExpression.toString(a)
      case terminal: Terminal => terminal.sourceString
    })

    lhs match {
      case Success(c:Call) =>
        c.task.outputs.find {_.name == rhs}.getOrElse {
          throw new SyntaxError(wdlSyntaxErrorFormatter.memberAccessReferencesBadTaskInput(ast))
        }
        Success(c)
      case f => f
    }
  }

  override def toRawString = ???

  validate()
}
