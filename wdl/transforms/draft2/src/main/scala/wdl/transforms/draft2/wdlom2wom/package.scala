package wdl.transforms.draft2

package object wdlom2wom {
  implicit val draft2WomCallableMaker: WdlDraft2WomCallableMaker.type =
    WdlDraft2WomCallableMaker
  implicit val draft2WomCallNodeMaker: WdlDraft2WomCallNodeMaker.type =
    WdlDraft2WomCallNodeMaker
  implicit val draft2WomConditionalNodeMaker: WdlDraft2WomConditionalNodeMaker.type =
    WdlDraft2WomConditionalNodeMaker
  implicit val draft2WomGraphMaker: WdlDraft2WomGraphMaker.type =
    WdlDraft2WomGraphMaker
  implicit val draft2WomScatterNodeMaker: WdlDraft2WomScatterNodeMaker.type =
    WdlDraft2WomScatterNodeMaker
  implicit val draft2WomTaskDefinitionMaker: WdlDraft2WomCommandTaskDefinitionMaker.type =
    WdlDraft2WomCommandTaskDefinitionMaker
  implicit val draft2WomWorkflowDefinitionMaker: WdlDraft2WomWorkflowDefinitionMaker.type =
    WdlDraft2WomWorkflowDefinitionMaker
}
