package wdl.transforms.draft2

package object wdlom2wom {
  implicit val draft2WomCallableMaker = WdlDraft2WomCallableMaker
  implicit val draft2WomCallNodeMaker = WdlDraft2WomCallNodeMaker
  implicit val draft2WomConditionalNodeMaker = WdlDraft2WomConditionalNodeMaker
  implicit val draft2WomGraphMaker = WdlDraft2WomGraphMaker
  implicit val draft2WomScatterNodeMaker = WdlDraft2WomScatterNodeMaker
  implicit val draft2WomTaskDefinitionMaker = WdlDraft2WomCommandTaskDefinitionMaker
  implicit val draft2WomWorkflowDefinitionMaker = WdlDraft2WomWorkflowDefinitionMaker
}
