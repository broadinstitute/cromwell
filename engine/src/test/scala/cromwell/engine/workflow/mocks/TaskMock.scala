package cromwell.engine.workflow.mocks

import cromwell.engine.workflow.mocks.DeclarationMock.DeclarationMockType
import org.specs2.mock.Mockito
import wdl.draft2.parser.WdlParser.Ast
import wdl.draft2.model.{Declaration, TaskOutput, WdlRuntimeAttributes, WdlTask}

trait TaskMock extends Mockito {
  
  def mockTask(name: String,
               declarations: Seq[Declaration] = Seq.empty,
               runtimeAttributes: WdlRuntimeAttributes = new WdlRuntimeAttributes(Map.empty),
               commandTemplateString: String = "!!shazam!!",
               outputs: Seq[DeclarationMockType] = Seq.empty
              ) = {
    val task = mock[WdlTask]
    task.declarations returns declarations
    task.runtimeAttributes returns runtimeAttributes
    task.commandTemplateString returns commandTemplateString
    task.name returns name
    task.unqualifiedName returns name
    task.outputs returns (outputs map {
      case (outputName, womType, expression) => TaskOutput(outputName, womType, expression, mock[Ast], Option(task))
    })
    task
  }
}
