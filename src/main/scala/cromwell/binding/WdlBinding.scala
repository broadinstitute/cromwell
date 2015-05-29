package cromwell.binding

import java.io.File

import cromwell.binding.WdlBinding.ImportResolver
import cromwell.binding.command._
import cromwell.binding.values.WdlValue
import cromwell.binding.types._
import cromwell.parser.WdlParser
import cromwell.parser.WdlParser._
import cromwell.util.FileUtil

import scala.collection.JavaConverters._
import scala.util.{Failure, Try}

/**
 * Main interface into the `cromwell.binding` package.
 *
 * Example usage:
 *
 * {{{
 * val binding = WdlBinding.process(new File("/path/to/file.wdl"))
 * binding.workflow.calls foreach { call =>
 *      println(call)
 * }
 * }}}
 */
object WdlBinding {

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
  def process(wdlFile: File): WdlBinding =
    process(readFile(wdlFile), wdlFile.toString, localImportResolver, None)

  def process(wdlFile: File, importResolver: ImportResolver): WdlBinding =
    process(readFile(wdlFile), wdlFile.toString, importResolver, None)

  def readFile(wdlFile: File): WdlSource = FileUtil.slurp(wdlFile)

  def process(wdlSource: WdlSource): WdlBinding =
    process(wdlSource, "string", localImportResolver, None)

  def process(wdlSource: WdlSource, importResolver: ImportResolver): WdlBinding =
    process(wdlSource, "string", importResolver, None)

  def process(wdlSource: WdlSource, resource: String): WdlBinding =
    process(wdlSource, resource, localImportResolver, None)

  private def process(wdlSource: WdlSource, resource: String, importResolver: ImportResolver, namespace: Option[String]): WdlBinding =
    new WdlBinding(WdlBinding.getAst(wdlSource, resource), wdlSource, importResolver, namespace)

  def getAst(wdlSource: WdlSource, resource: String): Ast = {
    val parser = new WdlParser()
    val tokens = parser.lex(wdlSource, resource)
    val terminalMap = (tokens.asScala.toVector map {(_, wdlSource)}).toMap
    val syntaxErrorFormatter = new WdlSyntaxErrorFormatter(terminalMap)
    parser.parse(tokens, syntaxErrorFormatter).toAst.asInstanceOf[Ast]
  }

  def process(wdlSource: WdlSource, resource: String, importResolver: ImportResolver): WdlBinding =
    process(wdlSource, resource, importResolver, None)

  /**
   * Given a WDL file, this will simply parse it and return the syntax tree
   * @param wdlFile The file to parse
   * @return an Abstract Syntax Tree (WdlParser.Ast) representing the structure of the code
   * @throws WdlParser.SyntaxError if there was a problem parsing the source code
   */
  def getAst(wdlFile: File): Ast = getAst(FileUtil.slurp(wdlFile), wdlFile.getName)

  def findAsts(ast: AstNode, name: String): Seq[Ast] = {
    ast match {
      case x: Ast =>
        val thisAst = if (x.getName.equals(name)) Seq(x) else Seq.empty[Ast]
        x.getAttributes.values.asScala.flatMap(findAsts(_, name)).toSeq ++ thisAst
      case x: AstList => x.asScala.toVector.flatMap(findAsts(_, name)).toSeq
      case x: Terminal => Seq.empty[Ast]
      case _ => Seq.empty[Ast]
    }
  }

  private def combine[T, U](map1: Map[T, Seq[U]], map2: Map[T, Seq[U]]): Map[T, Seq[U]] = {
    map1 ++ map2.map{ case (k,v) => k -> (v ++ map1.getOrElse(k, Seq.empty)) }
  }

  def findAstsWithTrail(ast: AstNode, name: String, trail: Seq[AstNode] = Seq.empty): Map[Ast, Seq[AstNode]] = {
    ast match {
      case x: Ast =>
        val thisAst = if (x.getName.equals(name)) Map(x -> trail) else Map.empty[Ast, Seq[AstNode]]
        combine(x.getAttributes.values.asScala.flatMap{y => findAstsWithTrail(y, name, trail :+ x)}.toMap, thisAst)
      case x: AstList => x.asScala.toVector.flatMap{y => findAstsWithTrail(y, name, trail :+ x)}.toMap
      case x: Terminal => Map.empty[Ast, Seq[AstNode]]
      case _ => Map.empty[Ast, Seq[AstNode]]
    }
  }

  def findTerminals(ast: AstNode): Seq[Terminal] = {
    ast match {
      case x: Ast => x.getAttributes.values.asScala.flatMap(findTerminals).toSeq
      case x: AstList => x.asScala.toVector.flatMap(findTerminals).toSeq
      case x: Terminal => Seq(x)
      case _ => Seq.empty[Terminal]
    }
  }

  def terminalMap(ast: Ast, source: WdlSource) = (findTerminals(ast) map {(_, source)}).toMap

  def getCallInput(ast: Ast): Map[String, WdlExpression] = getCallInputAsts(ast).map(processCallInput).toMap

  def getCallInputAsts(ast: Ast): Seq[Ast] = {
    findAsts(ast, AstNodeName.Inputs) match {
      case x: Seq[Ast] if x.size == 1 => WdlBinding.findAsts(x.head.getAttribute("map"), AstNodeName.IOMapping)
      case _ => Seq.empty[Ast]
    }
  }

  private def getTaskOutput(ast: Ast): TaskOutput = {
    val wdlType = getWdlType(ast.getAttribute("type"))
    val name = ast.getAttribute("var").asInstanceOf[Terminal].getSourceString
    val expression = ast.getAttribute("expression")
    new TaskOutput(name, wdlType, new WdlExpression(expression))
  }

  private def getCommand(astList: AstList): Command = {
    val parts = astList.asScala.toVector.map {
      case x: Terminal => new StringCommandPart(x.getSourceString)
      case x: Ast => getCommandParameter(x)
    }
    new Command(parts)
  }

  private def getCommandParameter(ast: Ast): ParameterCommandPart = {
    val wdlType = getWdlType(ast.getAttribute("type"))
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    new ParameterCommandPart(wdlType, name)
  }

  private def getWdlType(ast: AstNode): WdlType = {
    ast match {
      case t: Terminal if t.getSourceString == WdlFileType.toWdlString => WdlFileType
      case t: Terminal if t.getSourceString == WdlStringType.toWdlString => WdlStringType
      case t: Terminal if t.getSourceString == WdlIntegerType.toWdlString => WdlIntegerType
      case t: Terminal if t.getSourceString == WdlFloatType.toWdlString => WdlFloatType
      case t: Terminal if t.getSourceString == WdlBooleanType.toWdlString => WdlBooleanType
      case t: Terminal if t.getSourceString == WdlObjectType.toWdlString => WdlObjectType
      case null => WdlStringType
      case _ => throw new UnsupportedOperationException("Implement this later for compound types")
    }
  }

  private def processCallInput(ast: Ast): (String, WdlExpression) = {
    val key = ast.getAttribute("key").asInstanceOf[Terminal].getSourceString
    val expression = new WdlExpression(ast.getAttribute("value"))
    (key, expression)
  }

  object AstNodeName {
    val Task = "Task"
    val Workflow = "Workflow"
    val Command = "RawCommand"
    val Output = "Output"
    val CommandParameter = "CommandParameter"
    val Call = "Call"
    val IOMapping = "IOMapping"
    val Inputs = "Inputs"
    val MemberAccess = "MemberAccess"
  }

}

case class WdlBinding(ast: Ast, source: WdlSource, importResolver: ImportResolver, namespace: Option[String]) {

  import WdlBinding.AstNodeName

  val invalidTask = Task("INVALID", Command(Seq.empty[CommandPart]), Seq.empty[TaskOutput])

  /**
   * All `import` statement strings at the top of the document
   */
  val imports = ast.getAttribute("imports").asInstanceOf[AstList].asScala.toVector.map { i =>
    val uri = sourceString(i.asInstanceOf[Ast].getAttribute("uri"))
    val namespaceAst = i.asInstanceOf[Ast].getAttribute("namespace")
    val namespace = if (namespaceAst == null) None else Option(sourceString(namespaceAst))
    Import(uri, namespace)
  }

  /* WdlBinding objects for each import statement */
  val importedBindings = for {
    i <- imports
    source = importResolver(i.uri)
    if source.length > 0
  } yield WdlBinding.process(source, i.uri, importResolver, i.namespace)

  /* Create a map of Terminal -> WdlBinding */
  val terminalMap = WdlBinding.terminalMap(ast, source)
  val combinedTerminalMap = ((importedBindings map {x => x.terminalMap}) ++ Seq(terminalMap)) reduce (_ ++ _)
  val wdlSyntaxErrorFormatter = new WdlSyntaxErrorFormatter(combinedTerminalMap)

  /**
   * All imported `task` definitions for `import` statements without a namespace (e.g. no `as` clause)
   * These tasks are considered to be in this current workspace
   */
  val importedTaskAsts: Seq[Ast] = importedBindings flatMap { b =>
    b.namespace match {
      case None => b.taskAsts
      case _ => Seq.empty[Ast]
    }
  }
  val importedTasks: Seq[Task] = importedBindings flatMap { b =>
    b.namespace match {
      case None => b.tasks
      case _ => Seq.empty[Task]
    }
  }

  /**
   * All `task` definitions defined in the WDL file (i.e. not imported)
   */
  val localTaskAsts = WdlBinding.findAsts(ast, AstNodeName.Task)
  val localTasks = localTaskAsts.map(processTask)

  /**
   * All `task` definitions, including local and imported ones
   */
  val taskAsts = localTaskAsts ++ importedTaskAsts
  val tasks = localTasks ++ importedTasks

  /**
   * All imported `Workflow`s
   */
  val importedWorkflows: Seq[Workflow] = importedBindings flatMap { b =>
    b.namespace match {
      case None => b.workflows
      case _ => Seq.empty[Workflow]
    }
  }

  /**
   * All `Workflow`s defined in the WDL file (i.e. not imported)
   */
  val localWorkflows = WdlBinding.findAsts(ast, AstNodeName.Workflow).map(processWorkflow)

  /**
   * All `Workflow` definitions, including local and imported ones
   */
  val workflows = importedWorkflows ++ localWorkflows

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
      importedBindings find {_.namespace == Some(parts(0))} flatMap {_.findTaskAst(parts(1))}
    } else {
      taskAsts.find{t => sourceString(t.getAttribute("name")) == name}
    }
  }

  def findTask(name: String): Option[Task] = {
    if (name.contains(".")) {
      val parts = name.split("\\.", 2)
      importedBindings find {_.namespace == Some(parts(0))} flatMap {_.findTask(parts(1))}
    } else {
      tasks.find(_.name == name)
    }
  }

  private def processWorkflow(ast: Ast): Workflow = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val calls = WdlBinding.findAsts(ast, AstNodeName.Call).map {
      processCall
    }
    new Workflow(name, calls)
  }

  private def processTask(ast: Ast): Task = {
    val name = ast.getAttribute("name").asInstanceOf[Terminal].getSourceString
    val commandAsts = WdlBinding.findAsts(ast, AstNodeName.Command)
    if (commandAsts.size != 1) throw new UnsupportedOperationException("Expecting only one Command AST")
    val command = WdlBinding.getCommand(commandAsts.head.getAttribute("parts").asInstanceOf[AstList])
    val outputs = WdlBinding.findAsts(ast, AstNodeName.Output).map(WdlBinding.getTaskOutput)
    new Task(name, command, outputs)
  }

  private def processCall(ast: Ast): Call = {
    val alias: Option[String] = ast.getAttribute("alias") match {
      case x: Terminal => Option(x.getSourceString)
      case _ => None
    }
    val taskName = sourceString(ast.getAttribute("task"))
    val task = findTask(taskName) getOrElse invalidTask
    val inputs = WdlBinding.getCallInput(ast)
    new Call(alias, task, inputs.toMap)
  }

  /**
   * Confirm all required inputs are present and attempt to coerce raw inputs to `WdlValue`s.
   * This can fail if required raw inputs are missing or if the values for a specified raw input
   * cannot be coerced to the target type of the input as specified in the binding.
   */
  def coerceRawInputs(rawInputs: WorkflowRawInputs): Try[WorkflowCoercedInputs] = {

    def coerceRawInput(fqn: FullyQualifiedName, wdlType: WdlType): Try[WdlValue] = fqn match {
      case _ if rawInputs.contains(fqn) => wdlType.coerceRawValue(rawInputs.get(fqn).get)
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

  private def sourceString(astNode: AstNode): String = astNode.asInstanceOf[Terminal].getSourceString

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
    val workflowAsts = WdlBinding.findAsts(ast, AstNodeName.Workflow)

    if (workflowAsts.size > 1) {
      throw new SyntaxError(wdlSyntaxErrorFormatter.tooManyWorkflows(workflowAsts.asJava))
    }

    tasks foreach {task =>
      val taskAstsWithSameName = taskAsts filter {a => sourceString(a.getAttribute("name")) == task.name}
      /* No two tasks can have the same name */
      if (taskAstsWithSameName.size > 1) {
        throw new SyntaxError(wdlSyntaxErrorFormatter.duplicateTask(taskAstsWithSameName))
      }
      /* A task can not have duplicated inputs */
      val commandLineInputs = WdlBinding.findAsts(taskAstsWithSameName.head, AstNodeName.CommandParameter)
      commandLineInputs foreach {input =>
        val inputName = sourceString(input.getAttribute("name"))
        val inputsWithSameName = commandLineInputs filter {i => sourceString(i.getAttribute("name")) == inputName}
        if (inputsWithSameName.size > 1) {
          throw new SyntaxError(wdlSyntaxErrorFormatter.taskHasDuplicatedInputs(taskAstsWithSameName.head, inputsWithSameName))
        }
      }
    }

    /* Ensure that no namespaces collide with task/workflow names */
    ast.getAttribute("imports").asInstanceOf[AstList].asScala.toVector.map {i =>
      val namespaceAst = i.asInstanceOf[Ast].getAttribute("namespace").asInstanceOf[Terminal]
      if (namespaceAst != null) {
        findTaskAst(sourceString(namespaceAst)) match {
          case Some(taskAst) =>
            throw new SyntaxError(wdlSyntaxErrorFormatter.taskAndNamespaceHaveSameName(taskAst, namespaceAst))
          case _ =>
        }
        workflowAsts foreach {workflowAst =>
          if (sourceString(namespaceAst) == sourceString(workflowAst.getAttribute("name"))) {
            throw new SyntaxError(wdlSyntaxErrorFormatter.workflowAndNamespaceHaveSameName(workflowAst, namespaceAst))
          }
        }
      }
    }

    workflowAsts foreach { workflowAst =>
      val workflow = localWorkflows.head
      WdlBinding.findAsts(workflowAst, AstNodeName.Call) foreach { callAst =>
        val taskName = sourceString(callAst.getAttribute("task"))
        val taskAst = findTaskAst(taskName) getOrElse {
          throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskName(callAst, taskName))
        }
        val task = findTask(taskName) getOrElse {
          throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskName(callAst, taskName))
        }
        val callName = if(callAst.getAttribute("alias") != null) sourceString(callAst.getAttribute("alias")) else task.name
        val call = workflow.calls.find{_.name == callName} getOrElse {
          throw new SyntaxError(s"Cannot find call with name $callName")
        }

        call.inputMappings foreach { inputKv =>
          task.inputs find { taskInput => taskInput._1 == inputKv._1 } getOrElse {
            val callInput = WdlBinding.getCallInputAsts(callAst) find { p =>
              sourceString(p.getAttribute("key")) == inputKv._1
            } getOrElse {
              throw new SyntaxError(s"Can't find call input: ${inputKv._1}")
            }
            throw new SyntaxError(wdlSyntaxErrorFormatter.callReferencesBadTaskInput(callInput, taskAst))
          }

          /* All MemberAccess ASTs that are not contained in other MemberAccess ASTs */
          val memberAccesses = WdlBinding.findAstsWithTrail(inputKv._2.ast, "MemberAccess") filter {
            case (k, v) => !(v map {case a:Ast => a.getName} contains "MemberAccess")
          }
          memberAccesses.keys map validateMemberAccess
        }
      }
    }
  }

  /* Partially evaluate MemberAccess ASTs to make sure they're not nonsense at compile time */
  def validateMemberAccess(ast: Ast): Any = {
    val rhs = sourceString(ast.getAttribute("rhs"))
    val lhs = ast.getAttribute("lhs") match {
      case a: Ast => validateMemberAccess(a)
      case terminal: Terminal =>
        val name = sourceString(terminal)
        val bindings = importedBindings filter {b => b.namespace == Some(name)}
        val taskList = tasks filter {task => task.name == name}
        val all = bindings ++ taskList
        if (all.size == 0) throw new SyntaxError(wdlSyntaxErrorFormatter.undefinedMemberAccess(ast))
        all.head
    }

    lhs match {
      case b:WdlBinding =>
        val all = b.executables find {
          case w:Workflow => w.name == rhs
          case t:Task => t.name == rhs
        }
        if (all.size == 0) throw new SyntaxError(wdlSyntaxErrorFormatter.undefinedMemberAccess(ast))
        all.head
      case t:Task =>
        val wdlValue = t.outputs find{p => p.name == rhs}
        wdlValue getOrElse {
          throw new SyntaxError(wdlSyntaxErrorFormatter.memberAccessReferencesBadTaskInput(ast))
        }
    }
  }

  validate()
}
