package wdl.draft3.transforms

package object ast2wdlom {

  implicit val astFromAstNode = AstFromAstNode
  implicit val draft3FileElementFromAst = Draft3FileElementFromAst
  implicit val draft3ImportElementFromAst = Draft3ImportElementFromAst
  implicit val draft3TaskDefinitionElementFromAst = Draft3TaskDefinitionElementFromAst
  implicit val draft3WorkflowDefinitionElementFromAst = Draft3WorkflowDefinitionElementFromAst

}
