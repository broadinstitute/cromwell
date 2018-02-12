package wdl.draft3.transforms

package object ast2wdlom {

  implicit val astFromAstNode = AstFromAstNode
  implicit val draft3FileElementFromAstNode = FromAtoB.viaX(AstFromAstNode, Draft3FileElementFromAst)
  implicit val draft3ImportElementFromAstNode = FromAtoB.viaX(AstFromAstNode, Draft3ImportElementFromAst)
  implicit val draft3TaskDefinitionElementFromAstNode = FromAtoB.viaX(AstFromAstNode, Draft3TaskDefinitionElementFromAst)
  implicit val draft3WorkflowDefinitionElementFromAstNode = FromAtoB.viaX(AstFromAstNode, Draft3WorkflowDefinitionElementFromAst)

}
