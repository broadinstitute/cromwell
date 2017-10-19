package wdl

import java.nio.file.Path

import better.files._
import wdl4s.parser.WdlParser
import wdl4s.parser.WdlParser._
import wdl.WdlExpression.{AstForExpressions, AstNodeForExpressions}
import wdl.expression.ValueEvaluator.InterpolationTagPattern
import wom.core._
import wom.types._
import wom.values._

import scala.collection.JavaConverters._
import scala.language.postfixOps

object AstTools {

  class InterpolatedTerminal(val rootTerminal: Terminal, innerTerminal: Terminal, columnOffset: Int) extends Terminal(
    innerTerminal.getId,
    innerTerminal.getTerminalStr,
    innerTerminal.getSourceString,
    innerTerminal.getResource,
    // We want the line of the rootTerminal because innerTerminal is created purely from the expression string
    // So its first line will be line 0, which is very unlikely to be line 0 in the WDL file
    rootTerminal.getLine,
    // Column offset tells us the position of innerTerminal in the rootTerminal string
    // Therefore rootTerminal.getColumn + columnOffset gives us the column position for innerTerminal in the file (+1 to be on the $)
    rootTerminal.getColumn + columnOffset + 1
  )

  implicit class EnhancedAstNode(val astNode: AstNode) extends AnyVal {
    def findAsts(name: String): Seq[Ast] = AstTools.findAsts(astNode, name)
    def findAstsWithTrail(name: String, trail: Seq[AstNode] = Seq.empty): Map[Ast, Seq[AstNode]] = {
      astNode match {
        case x: Ast =>
          val thisAst = if (x.getName.equals(name)) Map(x -> trail) else Map.empty[Ast, Seq[AstNode]]
          combine(x.getAttributes.values.asScala.flatMap{_.findAstsWithTrail(name, trail :+ x)}.toMap, thisAst)
        case x: AstList => x.asScala.toVector.flatMap{_.findAstsWithTrail(name, trail :+ x)}.toMap
        case _: Terminal => Map.empty[Ast, Seq[AstNode]]
        case _ => Map.empty[Ast, Seq[AstNode]]
      }
    }

    private def findTerminalsInInterpolatedString(t: Terminal,
                                                  terminalType: String,
                                                  trail: Seq[AstNode],
                                                  parentTerminal: Option[Terminal],
                                                  columnOffset: Int) = {
      /*
        * Find all interpolations in the string terminal.
        * e.g: String a = "hello ${you}"
        * We'll create an expression from "you" and remember the position in the string 
        * "hello ${you}" at which we found "${you}".
       */
      val interpolatedExpressionAstNodesAndTheirMatchPosition = InterpolationTagPattern
        .findAllMatchIn(t.getSourceString)
        .foldLeft(List.empty[(AstNode, Int)])((nodes, exprValue) => {
          // This is the interpolated expression e.g ${my_var}
          val v = exprValue.group(0)
          // Create an expression from the content and remember the position of the match in the overall terminal string
          // so we can point to it in the error message if needed
          (WdlExpression.fromString(v.substring(2, v.length - 1)).ast, exprValue.start) +: nodes
        })

      interpolatedExpressionAstNodesAndTheirMatchPosition match {
        // If there's no interpolated expression and the parent terminal is of the right type,
        // create an interpolated terminal and we're done
        case Nil if t.getTerminalStr == terminalType =>
          val finalTerminal = parentTerminal map { parent => new InterpolatedTerminal(parent, t, columnOffset) } getOrElse t
          Map(finalTerminal -> trail)
        // No interpolated terminal and the parent terminal is not of the right type, we're done
        case Nil => Map.empty[Terminal, Seq[AstNode]]
        // We found some interpolated terminals, recursively find their inner terminals and propagate the root terminal.
        // Also propagate the accumulated columnOffset. The regex index match will start
        // over at 0 in the next round of matching so we need to keep track of the offset as we recurse
        case expressions => expressions.flatMap({
          case (innerNode, offset) => innerNode.findTerminalsWithTrail(terminalType, trail :+ t, Option(parentTerminal.getOrElse(t)), columnOffset + offset)
        }).toMap
      }
    }

    def findTerminalsWithTrail(terminalType: String,
                               trail: Seq[AstNode] = Seq.empty,
                               parentTerminal: Option[Terminal] = None,
                               columnOffset: Int = 0): Map[Terminal, Seq[AstNode]] = {
      astNode match {
        case a: Ast => a.getAttributes.values.asScala flatMap { _.findTerminalsWithTrail(terminalType, trail :+ a) } toMap
        case a: AstList => a.asScala.toVector flatMap { _.findTerminalsWithTrail(terminalType, trail :+ a) } toMap
        case t: Terminal if t.getTerminalStr == terminalType =>
          val finalTerminal = parentTerminal map { parent => new InterpolatedTerminal(parent, t, columnOffset) } getOrElse t
          Map(finalTerminal -> trail)
        case t: Terminal => findTerminalsInInterpolatedString(t, terminalType, trail, parentTerminal, columnOffset)
        case _ => Map.empty[Terminal, Seq[AstNode]]
      }
    }
    def findFirstTerminal: Option[Terminal] = {
      Option(astNode) flatMap {
        case l: AstList => l.astListAsVector.flatMap(_.findFirstTerminal).headOption
        case a: Ast => a.getAttributes.asScala.toMap.flatMap({ case (_, v) => v.findFirstTerminal }).headOption
        case t: Terminal => Option(t)
      }
    }
    def findTopLevelMemberAccesses(): Iterable[Ast] = AstTools.findTopLevelMemberAccesses(astNode)
    def sourceString: String = astNode.asInstanceOf[Terminal].getSourceString
    def astListAsVector: Seq[AstNode] = astNode.asInstanceOf[AstList].asScala.toVector
    def wdlType(wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): WdlType = {
      astNode match {
        case t: Terminal =>
          t.getSourceString match {
            case WdlFileType.toWdlString => WdlFileType
            case WdlStringType.toWdlString => WdlStringType
            case WdlIntegerType.toWdlString => WdlIntegerType
            case WdlFloatType.toWdlString => WdlFloatType
            case WdlBooleanType.toWdlString => WdlBooleanType
            case WdlObjectType.toWdlString => WdlObjectType
            case "Array" => throw new SyntaxError(wdlSyntaxErrorFormatter.arrayMustHaveATypeParameter(t))
          }
        case a: Ast if isOptionalType(a) => optionalType(a, wdlSyntaxErrorFormatter)
        case a: Ast if isNonEmptyType(a) => nonEmptyType(a, wdlSyntaxErrorFormatter)
        case a: Ast =>
          val subtypes = a.getAttribute("subtype").astListAsVector
          val typeTerminal = a.getAttribute("name").asInstanceOf[Terminal]
          a.getAttribute("name").sourceString match {
            case "Pair" =>
              if (subtypes.size != 2) throw new SyntaxError(wdlSyntaxErrorFormatter.pairMustHaveExactlyTwoTypeParameters(typeTerminal))
              val leftType = subtypes.head.wdlType(wdlSyntaxErrorFormatter)
              val rightType = subtypes.tail.head.wdlType(wdlSyntaxErrorFormatter)
              WdlPairType(leftType, rightType)
            case "Array" =>
              if (subtypes.size != 1) throw new SyntaxError(wdlSyntaxErrorFormatter.arrayMustHaveOnlyOneTypeParameter(typeTerminal))
              val member = subtypes.head.wdlType(wdlSyntaxErrorFormatter)
              WdlArrayType(member)
            case "Map" =>
              if (subtypes.size != 2) throw new SyntaxError(wdlSyntaxErrorFormatter.mapMustHaveExactlyTwoTypeParameters(typeTerminal))
              val keyType = subtypes.head.wdlType(wdlSyntaxErrorFormatter)
              val valueType = subtypes.tail.head.wdlType(wdlSyntaxErrorFormatter)
              WdlMapType(keyType, valueType)
          }
        case _ => throw new UnsupportedOperationException(s"Unexpected WDL type AST: ${astNode.sourceString}")
      }
    }

    def isOptionalType(a: Ast) = a.getName.equals("OptionalType")
    def optionalType(a: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter) = WdlOptionalType(a.getAttribute("innerType").wdlType(wdlSyntaxErrorFormatter))

    def isNonEmptyType(a: Ast) = a.getName.equals("NonEmptyType")
    def nonEmptyType(a: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter) = {
      val innerType = a.getAttribute("innerType").wdlType(wdlSyntaxErrorFormatter)
      innerType match {
        case arrayType: WdlArrayType => arrayType.asNonEmptyArrayType
        case _ => throw new UnsupportedOperationException("Currently the only supported non-empty types are Array[X]+")
      }
    }

    def wdlValue(wdlType: WdlType, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): WdlValue = {

      def astToMap(ast: Ast) = {
        val mapType = wdlType.asInstanceOf[WdlMapType]
        val elements = ast.getAttribute("map").asInstanceOf[AstList].asScala.toVector.map({ kvnode =>
          val k = kvnode.asInstanceOf[Ast].getAttribute("key").wdlValue(mapType.keyType, wdlSyntaxErrorFormatter)
          val v = kvnode.asInstanceOf[Ast].getAttribute("value").wdlValue(mapType.valueType, wdlSyntaxErrorFormatter)
          k -> v
        }).toMap

        WdlMap(mapType, elements)
      }

      def astToObject(ast: Ast) = {
        val elements = ast.getAttribute("map").asInstanceOf[AstList].asScala.toVector.map({ kvnode =>
          val k = kvnode.asInstanceOf[Ast].getAttribute("key").sourceString
          val v = kvnode.asInstanceOf[Ast].getAttribute("value").wdlValue(WdlStringType, wdlSyntaxErrorFormatter)
          k -> v
        }).toMap

        WdlObject(elements)
      }

      def astTupleToValue(a: Ast): WdlValue = {
        val subElements = a.getAttribute("values").astListAsVector
        if (subElements.size == 1) {
          // Tuple 1 is equivalent to the value inside it. Enables nesting parens, e.g. (1 + 2) + 3
          a.wdlValue(wdlType, wdlSyntaxErrorFormatter)
        } else if (subElements.size == 2 && wdlType.isInstanceOf[WdlPairType]) {
          val pairType = wdlType.asInstanceOf[WdlPairType]
          WdlPair(subElements.head.wdlValue(pairType.leftType, wdlSyntaxErrorFormatter), subElements(1).wdlValue(pairType.rightType, wdlSyntaxErrorFormatter))
        } else {
          throw new SyntaxError(s"Could not convert AST to a $wdlType (${Option(astNode).getOrElse("No AST").toString})")
        }
      }

      astNode match {
        case t: Terminal if t.getTerminalStr == "string" && wdlType == WdlStringType => WdlString(t.getSourceString)
        case t: Terminal if t.getTerminalStr == "string" && wdlType == WdlFileType => WdlFile(t.getSourceString)
        case t: Terminal if t.getTerminalStr == "integer" && wdlType == WdlIntegerType => WdlInteger(t.getSourceString.toInt)
        case t: Terminal if t.getTerminalStr == "float" && wdlType == WdlFloatType => WdlFloat(t.getSourceString.toDouble)
        case t: Terminal if t.getTerminalStr == "boolean" && wdlType == WdlBooleanType => t.getSourceString.toLowerCase match {
          case "true" => WdlBoolean.True
          case "false" => WdlBoolean.False
        }
        // TODO: The below cases, ArrayLiteral and MapLiteral, ObjectLiteral are brittle. They recursively call this wdlValue().
        // However, those recursive calls might contain full-on expressions instead of just other literals.  This
        // whole thing ought to be part of the regular expression evaluator, though I imagine that's non-trivial.
        case a: Ast if a.getName == "ArrayLiteral" && wdlType.isInstanceOf[WdlArrayType] =>
          val arrType = wdlType.asInstanceOf[WdlArrayType]
          val elements = a.getAttribute("values").astListAsVector map {node => node.wdlValue(arrType.memberType, wdlSyntaxErrorFormatter)}
          WdlArray(arrType, elements)
        case a: Ast if a.getName == "TupleLiteral" => astTupleToValue(a)
        case a: Ast if a.getName == "MapLiteral" && wdlType.isInstanceOf[WdlMapType] => astToMap(a)
        case a: Ast if a.getName == "ObjectLiteral" && wdlType == WdlObjectType => astToObject(a)
        case _ => throw new SyntaxError(s"Could not convert AST to a $wdlType (${Option(astNode).getOrElse("No AST").toString})")
      }
    }
  }

  implicit class EnhancedAstSeq(val astSeq: Seq[Ast]) extends AnyVal {
    def duplicatesByName: Seq[Ast] = {
      astSeq.groupBy(_.getAttribute("name").sourceString).collect({case (_ ,v) if v.size > 1 => v.head}).toVector
    }
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
    val Runtime = "Runtime"
    val RuntimeAttribute = "RuntimeAttribute"
    val Declaration = "Declaration"
    val WorkflowOutputWildcard = "WorkflowOutputWildcard"
    val WorkflowOutputDeclaration = "WorkflowOutputDeclaration"
    val WorkflowOutputs = "WorkflowOutputs"
    val Scatter = "Scatter"
    val Meta = "Meta"
    val ParameterMeta = "ParameterMeta"
    val Namespace = "Namespace" // TODO: rename this in the grammar
    val If = "If"
  }

  def getAst(workflowSource: WorkflowSource, resource: String): Ast = {
    val parser = new WdlParser()
    val tokens = parser.lex(workflowSource, resource)
    val terminalMap = (tokens.asScala.toVector map {(_, workflowSource)}).toMap
    val syntaxErrorFormatter = WdlSyntaxErrorFormatter(terminalMap)
    parser.parse(tokens, syntaxErrorFormatter).toAst.asInstanceOf[Ast]
  }

  /**
    * Given a WDL file, this will simply parse it and return the syntax tree
    *
    * @param wdlFile The file to parse
    * @return an Abstract Syntax Tree (WdlParser.Ast) representing the structure of the code
    * @throws WdlParser.SyntaxError if there was a problem parsing the source code
    */
  def getAst(wdlFile: Path): Ast = getAst(File(wdlFile).contentAsString, File(wdlFile).name)

  def findAsts(ast: AstNode, name: String): Seq[Ast] = {
    ast match {
      case x: Ast =>
        val thisAst = if (x.getName.equals(name)) Seq(x) else Seq.empty[Ast]
        x.getAttributes.values.asScala.flatMap(findAsts(_, name)).toSeq ++ thisAst
      case x: AstList => x.asScala.toVector.flatMap(findAsts(_, name))
      case _: Terminal => Seq.empty[Ast]
      case _ => Seq.empty[Ast]
    }
  }

  def findTerminals(ast: AstNode): Seq[Terminal] = {
    ast match {
      case x: Ast => x.getAttributes.values.asScala.flatMap(findTerminals).toSeq
      case x: AstList => x.asScala.toVector.flatMap(findTerminals)
      case x: Terminal => Seq(x)
      case _ => Seq.empty[Terminal]
    }
  }

  /**
    * All MemberAccess ASTs that are not contained in other MemberAccess ASTs
    *
    * The reason this returns a collection would be expressions such as "a.b.c + a.b.d", each one of those
    * would have its own MemberAccess - "a.b.c" and "a.b.d"
    */
  def findTopLevelMemberAccesses(expr: AstNode): Iterable[Ast] = expr.findAstsWithTrail("MemberAccess").filterNot {
    case (_, v) => v exists {
      case a: Ast => a.getName == "MemberAccess"
      case _ => false
    }
  }.keys

  final case class VariableReference(terminal: Terminal, trail: Iterable[AstNode]) {
    /**
      * If this is a simple MemberAccess (both sides are terminals),
      * find the rhs corresponding to "terminal" in the trail.
      * e.g
      *
      * if terminal is "a" and the trail contains a MemberAccess(lhs: Terminal("a"), rhs: Terminal("b")),
      * this will return Some(Terminal("b"))
      */
    lazy val terminalSubIdentifier: Option[Terminal] = trail.collectFirst {
      case a: Ast if a.isMemberAccess
        && a.getAttribute("lhs") == terminal
        && a.getAttribute("rhs").isTerminal => a.getAttribute("rhs").asInstanceOf[Terminal]
    }

    def fullVariableReferenceString = terminal.getSourceString + (terminalSubIdentifier map { "." + _.getSourceString } getOrElse "")
  }

  /**
    * All variable references in the expression AstNode that are not part of MemberAccess ASTs
    *
    * These represent anything that would need to be have scope resolution done on it to determine the value
    */
  def findVariableReferences(expr: AstNode): Iterable[VariableReference] = {
    def isMemberAccessRhs(identifier: Terminal, trail: Seq[AstNode]): Boolean = {
      // e.g. for MemberAccess ast representing source code A.B.C, this would return true for only B,C and not A
      trail.collect({ case a: Ast if a.isMemberAccess && a.getAttribute("rhs") == identifier => a }).nonEmpty
    }
    def isFunctionName(identifier: Terminal, trail: Seq[AstNode]): Boolean = {
      trail.lastOption match {
        case Some(last: Ast) if last.isFunctionCall && last.getAttribute("name") == identifier => true
        case _ => false
      }
    }

    /* terminal is the "lefter" lhs
     * trail is how we arrived to identifier from the original ast
     * e.g #1 (in "pseudo ast code"):
     * 
     * If MemberAccess is "a.b"
     * terminal will be Terminal("a")
     * trail will be Seq(
     *   MemberAccess(
     *     lhs: Terminal("a"),
     *     rhs: Terminal("b")
     *   )  
     * )
     * 
     * e.g #2:
     * If MemberAccess is "a.b.c"
     * terminal will be Terminal("a")
     * trail will be Seq(
     *   MemberAccess(
     *     lhs: MemberAccess(lhs: Terminal("a"), rhs: Terminal("b")),
     *     rhs: Terminal("c")
     *   ),
     *   MemberAccess(
     *     lhs: Terminal("a"),
     *     rhs: Terminal("b")
     *   )
     * )
     * 
     * There also might be other types of nodes in trail than MemberAccess depending the expression.
     */
    expr.findTerminalsWithTrail("identifier").collect({
      case (terminal, trail) if !isMemberAccessRhs(terminal, trail) && !isFunctionName(terminal, trail) => VariableReference(terminal, trail)
    })
  }

  /**
    * Given a Call AST, this will validate that there is 0 or 1 'input' section with non-empty
    * key/value pairs (if the input section is specified).  This will then return a Seq of
    * IOMapping(key=<terminal> value=<expression ast>)
    *
    * @param ast The call AST
    * @param wdlSyntaxErrorFormatter The wdl syntax error formatter
    * @return Seq[Ast] where the AST is a IOMapping(key=<terminal> value=<expression ast>)
    */
  def callInputSectionIOMappings(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Seq[Ast] = {
    val callTaskName = ast.getAttribute("task").asInstanceOf[Terminal]

    /* Filter out all empty 'input' sections first */
    val callInputSections = ast.findAsts(AstNodeName.Inputs).map {inputSectionAst =>
      inputSectionAst.getAttribute("map").findAsts(AstNodeName.IOMapping) match {
        case kvPairs: Seq[Ast] if kvPairs.isEmpty => throw new SyntaxError(wdlSyntaxErrorFormatter.emptyInputSection(callTaskName))
        case _ => inputSectionAst
      }
    }

    /* Then, make sure there is at most one 'input' section defined, then return the a Seq of IOMapping(key=<terminal>, value=<expression>) ASTs*/
    callInputSections match {
      case asts: Seq[Ast] if asts.size == 1 => asts.head.getAttribute("map").findAsts(AstNodeName.IOMapping)
      case asts: Seq[Ast] if asts.isEmpty => Seq.empty[Ast]
      case asts: Seq[Ast] =>
        /* Uses of .head here are assumed by the above code that ensures that there are no empty maps */
        val secondInputSectionIOMappings = asts(1).getAttribute("map").astListAsVector
        val firstKeyTerminal = secondInputSectionIOMappings.head.asInstanceOf[Ast].getAttribute("key").asInstanceOf[Terminal]
        throw new SyntaxError(wdlSyntaxErrorFormatter.multipleInputStatementsOnCall(firstKeyTerminal))
    }
  }

  def terminalMap(ast: Ast, source: WorkflowSource) = (findTerminals(ast) map {(_, source)}).toMap

  def wdlSectionToStringMap(ast: Ast, node: String, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): Map[String, String] = {
    ast.findAsts(node) match {
      case a if a.isEmpty => Map.empty[String, String]
      case a if a.size == 1 =>
        // Yes, even 'meta {}' and 'parameter_meta {}' sections have RuntimeAttribute ASTs.
        // In hindsight, this was a poor name for the AST.
        a.head.findAsts(AstNodeName.RuntimeAttribute).map({ ast =>
          val key = ast.getAttribute("key").asInstanceOf[Terminal]
          val value = ast.getAttribute("value")
          if (!value.isInstanceOf[Terminal] || value.asInstanceOf[Terminal].getTerminalStr != "string") {
            // Keys are parsed as identifiers, but values are parsed as expressions.
            // For now, only accept expressions that are strings
            throw new SyntaxError(wdlSyntaxErrorFormatter.expressionExpectedToBeString(key))
          }
          key.sourceString -> value.sourceString
        }).toMap
      case _ => throw new SyntaxError(wdlSyntaxErrorFormatter.expectedAtMostOneSectionPerTask(node, ast.getAttribute("name").asInstanceOf[Terminal]))
    }
  }

  private def combine[T, U](map1: Map[T, Seq[U]], map2: Map[T, Seq[U]]): Map[T, Seq[U]] = {
    map1 ++ map2.map{ case (k,v) => k -> (v ++ map1.getOrElse(k, Seq.empty)) }
  }
}
