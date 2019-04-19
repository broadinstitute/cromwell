package wdl.transforms.base.ast2wdlom

import common.Checked
import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.model.draft3.elements._
import wom.types._
import wdl.transforms.base.wdlom2wdl.WdlWriter.ops._
import wdl.transforms.base.wdlom2wdl.WdlWriterImpl.typeElementWriter

object AstNodeToTypeElement {

  def astNodeToTypeElement(additionalPrimitiveTypes: Map[String, WomPrimitiveType]): CheckedAtoB[GenericAstNode, TypeElement] =
    CheckedAtoB.fromCheck("convert AstNode to TypeElement") { astNode =>

      implicit lazy val astNodeToTypeElementInst = astNodeToTypeElement(additionalPrimitiveTypes)

      lazy val fullTypeMap = typeMap ++ additionalPrimitiveTypes

      astNode match {
        case a: GenericAst if a.getName == "OptionalType" => a.getAttributeAs[TypeElement]("innerType") map OptionalTypeElement
        case a: GenericAst if a.getName == "NonEmptyType" => a.getAttributeAs[TypeElement]("innerType") map NonEmptyTypeElement
        case a: GenericAst if a.getName == "Type" => compoundType(a)
        case unknownAst: GenericAst => s"No rule available to create TypeElement from Ast: '${unknownAst.getName}'".invalidNelCheck
        case t: GenericTerminal if t.getTerminalStr == "type" && fullTypeMap.contains(t.getSourceString) => PrimitiveTypeElement(fullTypeMap(t.getSourceString)).validNelCheck
        case t: GenericTerminal if t.getTerminalStr == "type" && t.getSourceString == "Object" => ObjectTypeElement.validNelCheck
        case t: GenericTerminal if t.getTerminalStr == "identifier" => TypeAliasElement(t.getSourceString).validNelCheck
        case t: GenericTerminal => s"No rule available to create TypeElement from '${t.getTerminalStr}' Terminal with value '${t.getSourceString}'".invalidNelCheck
        case _ => s"No rule available to create TypeElement from AstNode: ${astNode.getClass.getSimpleName}".invalidNelCheck
      }
    }

  private def compoundType(typeAst: GenericAst)
                          (implicit astNodeToExpressionElement: CheckedAtoB[GenericAstNode, TypeElement]): Checked[TypeElement] = typeAst.getAttributeAs[String]("name") flatMap {
    case "Array" => typeAst.getAttributeAsVector[TypeElement]("subtype") flatMap {
      case one if one.size == 1 => ArrayTypeElement(one.head).validNelCheck
      case other => s"Arrays must have exactly one type parameter, but got ${other.map(_.toWdlV1).mkString("[", ",", "]")}".invalidNelCheck
    }
    case "Pair" => typeAst.getAttributeAsVector[TypeElement]("subtype") flatMap {
      case two if two.size == 2 => PairTypeElement(two.head, two(1)).validNelCheck
      case other => s"Pairs must have exactly two type parameters, but got ${other.map(_.toWdlV1).mkString("[", ",", "]")}".invalidNelCheck
    }
    case "Map" => typeAst.getAttributeAsVector[TypeElement]("subtype") flatMap {
      case two if two.size == 2 => MapTypeElement(two.head, two(1)).validNelCheck
      case other => s"Maps must have exactly two type parameters, but got ${other.map(_.toWdlV1).mkString("[", ",", "]")}".invalidNelCheck
    }
    case unknown => s"No rule available to create TypeElement from compound type: $unknown".invalidNelCheck
  }

  private val typeMap: Map[String, WomPrimitiveType] = Map(
    "Int" -> WomIntegerType,
    "String" -> WomStringType,
    "Float" -> WomFloatType,
    "Boolean" -> WomBooleanType,
    "File" -> WomSingleFileType
  )
}
