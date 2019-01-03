package wdl.model.draft3.elements

/**
  * A supertype for elements in the language which represent nodes in the workflow AST.
  */
trait WorkflowGraphElement extends LanguageElement with WorkflowBodyElement
