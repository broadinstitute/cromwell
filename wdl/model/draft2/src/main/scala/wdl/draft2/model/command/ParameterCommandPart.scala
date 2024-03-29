package wdl.draft2.model.command

import cats.syntax.validated._
import common.validation.ErrorOr._
import wdl.draft2.model.AstTools.EnhancedAstNode
import wdl.draft2.model.exception.VariableNotFoundException
import wdl.draft2.model.expression.WdlFunctions
import wdl.draft2.model.{Declaration, WdlExpression, WdlSyntaxErrorFormatter}
import wdl.draft2.parser.WdlParser.{Ast, SyntaxError, Terminal}
import wom.types.WomOptionalType
import wom.values._
import wom.{CommandSetupSideEffectFile, InstantiatedCommand, OptionalNotSuppliedException}

import scala.language.postfixOps
import scala.util.{Failure, Success}

object ParameterCommandPart {
  def apply(ast: Ast, wdlSyntaxErrorFormatter: WdlSyntaxErrorFormatter): ParameterCommandPart = {
    val attributes = ast.getAttribute("attributes").astListAsVector map { a =>
      val ast = a.asInstanceOf[Ast]
      (ast.getAttribute("key").sourceString, ast.getAttribute("value").sourceString)
    } toMap
    val expression = WdlExpression(ast.getAttribute("expr"))
    if (
      (attributes.contains("true") && !attributes.contains("false")) || (attributes.contains("false") && !attributes
        .contains("true"))
    ) {
      // .head because we can't get here without there being at least one attribute
      val firstAttr =
        ast.getAttribute("attributes").astListAsVector.head.asInstanceOf[Ast].getAttribute("key").asInstanceOf[Terminal]
      throw new SyntaxError(wdlSyntaxErrorFormatter.trueAndFalseAttributesAreRequired(firstAttr))
    }
    new ParameterCommandPart(attributes, expression)
  }
}

case class ParameterCommandPart(attributes: Map[String, String], expression: WdlExpression) extends WdlCommandPart {
  def attributesToString: String = if (attributes.nonEmpty)
    attributes.map { case (k, v) => s"$k=${WomString(v).toWomString}" }.mkString(" ") + " "
  else ""
  override def toString: String = "${" + s"$attributesToString${expression.toWomString}" + "}"

  override def instantiate(declarations: Seq[Declaration],
                           inputs: Map[String, WomValue],
                           functions: WdlFunctions[WomValue],
                           valueMapper: (WomValue) => WomValue
  ): ErrorOr[List[InstantiatedCommand]] = {
    // This is a safety net.
    // In Cromwell's production code, optional declarations are always passed to instantiate, as WdlOptionalValue.none(type) if necessary.
    def lookupDeclaration(s: String) = declarations.collectFirst {
      // The backtick syntax (`s`) allows us to equality-check 's' against the match/case result:
      case Declaration(womType: WomOptionalType, `s`, _, _, _) => womType.none
    } getOrElse { throw VariableNotFoundException(s) }

    val lookup: String => WomValue = (s: String) => valueMapper(inputs.getOrElse(s, lookupDeclaration(s)))

    // Note that evaluating this expression may have the side effect of causing a file to be created if `writeFile`
    // is invoked. That file will be written to an "engine-relative" location which may be in cloud storage.
    val evaluatedCommandPartExpression: ErrorOr[WomValue] = expression.evaluate(lookup, functions) match {
      case Success(v) =>
        v match {
          case WomOptionalValue(_, opt) => opt.getOrElse(defaultString).validNel
          case _ => v.validNel
        }
      case Failure(OptionalNotSuppliedException(_)) => defaultString.validNel
      case Failure(e) => s"Could not evaluate expression: ${expression.toWomString}: ${e.getMessage}".invalidNel
    }

    // Create the stringified version of the command and record any file created in the process.
    def instantiateCommand(value: WomValue): ErrorOr[InstantiatedCommand] =
      (valueMapper(value), value) match {
        case (b: WomBoolean, _) if attributes.contains("true") && attributes.contains("false") =>
          InstantiatedCommand(if (b.value) attributes.get("true").head else attributes.get("false").head).validNel
        case (f: WomFile, unmappedFile: WomFile) if unmappedFile.value != f.value =>
          // Files generated by writeFiles have "engine-relative" paths which will be different from the container paths
          // calculated by `valueMapper`. "engine-relative" may mean either the non-Docker container path on the host
          // running Cromwell, or a cloud path. Capture these newly created files and their engine paths here.
          InstantiatedCommand(commandString = f.valueString,
                              createdFiles = List(CommandSetupSideEffectFile(unmappedFile))
          ).validNel
        case (f: WomFile, _) =>
          InstantiatedCommand(f.valueString).validNel
        case (p: WomPrimitive, _) => InstantiatedCommand(p.valueString).validNel
        case (a: WomArray, _) if attributes.contains("sep") =>
          InstantiatedCommand(a.value.map(_.valueString).mkString(attributes.get("sep").head)).validNel
        case (_: WomArray, _) =>
          s"Expression '${expression.toWomString}' evaluated to an Array but no 'sep' was specified".invalidNel
        case _ =>
          s"Could not string-ify value: $value".invalidNel
      }

    for {
      value <- evaluatedCommandPartExpression
      instantiatedCommand <- instantiateCommand(value)
    } yield List(instantiatedCommand)
  }

  lazy val defaultString = WomString(attributes.getOrElse("default", ""))
}
