package cromwell.parser

import java.io.File

import cromwell.binding._
import cromwell.binding.types._
import cromwell.parser.WdlParser._
import cromwell.util.FileUtil

import scala.collection.JavaConverters._

object AstTools {
  implicit class EnhancedAstNode(val astNode: AstNode) extends AnyVal {
    def findAsts(name: String): Seq[Ast] = AstTools.findAsts(astNode, name)
    def findAstsWithTrail(name: String, trail: Seq[AstNode] = Seq.empty): Map[Ast, Seq[AstNode]] = {
      AstTools.findAstsWithTrail(astNode, name, trail)
    }
    def findTerminals(): Seq[Terminal] = AstTools.findTerminals(astNode)
    def findTopLevelMemberAccesses(): Iterable[Ast] = AstTools.findTopLevelMemberAccesses(astNode)
    def sourceString(): String = astNode.asInstanceOf[Terminal].getSourceString
    
    def wdlType: WdlType = {
      astNode match {
        case t: Terminal =>
          t.getSourceString match {
            case WdlFileType.toWdlString => WdlFileType
            case WdlStringType.toWdlString => WdlStringType
            case WdlIntegerType.toWdlString => WdlIntegerType
            case WdlFloatType.toWdlString => WdlFloatType
            case WdlBooleanType.toWdlString => WdlBooleanType
            case WdlObjectType.toWdlString => WdlObjectType
          }
        case a: Ast =>
          val subtypes = a.getAttribute("subtype").asInstanceOf[AstList].asScala.toSeq
          val typeTerminal = a.getAttribute("name").asInstanceOf[Terminal]
          a.getAttribute("name").sourceString() match {
            case "Array" =>
              val member = subtypes.head.wdlType
              WdlArrayType(member)
          }
        case null => WdlStringType
        case _ => throw new UnsupportedOperationException("Implement this later for compound types")
      }
    }
  }

  implicit class EnhancedAstSeq(val astSeq: Seq[Ast]) extends AnyVal {
    def duplicatesByName: Seq[Ast] = {
      astSeq.groupBy(_.getAttribute("name").sourceString()).collect({case (_ ,v) if v.size > 1 => v.head}).toVector
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
    val Declaration = "Declaration"
  }

  def getAst(wdlSource: WdlSource, resource: String): Ast = {
    val parser = new WdlParser()
    val tokens = parser.lex(wdlSource, resource)
    val terminalMap = (tokens.asScala.toVector map {(_, wdlSource)}).toMap
    val syntaxErrorFormatter = new WdlSyntaxErrorFormatter(terminalMap)
    parser.parse(tokens, syntaxErrorFormatter).toAst.asInstanceOf[Ast]
  }

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

  def findAstsWithTrail(ast: AstNode, name: String, trail: Seq[AstNode] = Seq.empty): Map[Ast, Seq[AstNode]] = {
    ast match {
      case x: Ast =>
        val thisAst = if (x.getName.equals(name)) Map(x -> trail) else Map.empty[Ast, Seq[AstNode]]
        combine(x.getAttributes.values.asScala.flatMap{_.findAstsWithTrail(name, trail :+ x)}.toMap, thisAst)
      case x: AstList => x.asScala.toVector.flatMap{_.findAstsWithTrail(name, trail :+ x)}.toMap
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

  /* All MemberAccess ASTs that are not contained in other MemberAccess ASTs */
  def findTopLevelMemberAccesses(expr: AstNode): Iterable[Ast] = expr.findAstsWithTrail("MemberAccess").filterNot {
    case(k, v) => v.exists{case a:Ast => a.getName == "MemberAccess"}
  }.keys

  def callInputAsts(ast: Ast): Seq[Ast] = {
    findAsts(ast, AstNodeName.Inputs) match {
      case x: Seq[Ast] if x.size == 1 => x.head.getAttribute("map").findAsts(AstNodeName.IOMapping)
      case _ => Seq.empty[Ast]
    }
  }

  def terminalMap(ast: Ast, source: WdlSource) = (findTerminals(ast) map {(_, source)}).toMap

  private def combine[T, U](map1: Map[T, Seq[U]], map2: Map[T, Seq[U]]): Map[T, Seq[U]] = {
    map1 ++ map2.map{ case (k,v) => k -> (v ++ map1.getOrElse(k, Seq.empty)) }
  }
}

