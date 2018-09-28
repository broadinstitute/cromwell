package wdl.transforms.base.ast2wdlom

import common.transforms.CheckedAtoB
import common.validation.Checked._
import wdl.model.draft3.elements._

object AstNodeToStaticString {
  def astNodeToStaticStringElement(): CheckedAtoB[GenericAstNode, StaticString] = CheckedAtoB.fromCheck("convert AstNode to StaticString") {
    case a: GenericAst if a.getName == "StaticString" =>
      // The 'value' field (representing content between the "") is optional in the grammar.
      // No content between the "" represents an empty string:
      val value = if (a.getAttributes.contains("value")) {
        a.getAttributeAs[String]("value")
      } else "".validNelCheck
      value map StaticString.apply
    case other => s"Bad value $other".invalidNelCheck
  }
}
