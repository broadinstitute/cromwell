package cromwell.engine.workflow.mocks

import cromwell.engine.workflow.mocks.DeclarationMock.DeclarationMockType
import org.specs2.mock.Mockito
import wdl4s._
import wdl4s.parser.WdlParser.Ast

trait TaskMock extends Mockito {
  
  def mockTask(name: String,
               declarations: Seq[Declaration] = Seq.empty,
               runtimeAttributes: RuntimeAttributes = new RuntimeAttributes(Map.empty),
               commandTemplateString: String = "!!shazam!!",
               outputs: Seq[DeclarationMockType] = Seq.empty
              ) = {
    val task = mock[Task]
    task.declarations returns declarations
    task.runtimeAttributes returns runtimeAttributes
    task.commandTemplateString returns commandTemplateString
    task.name returns name
    task.unqualifiedName returns name
    task.outputs returns (outputs map {
      case (outputName, wdlType, expression) => TaskOutput(outputName, wdlType, expression, mock[Ast], Option(task))
    })
    task
  }
}
