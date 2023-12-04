package wdl.transforms.base.wdlom2wom.graph

import wdl.transforms.base.wdlom2wom.expression.renaming.IdentifierLookupRenamer
import wdl.model.draft3.elements._
import wdl.transforms.base.wdlom2wom.graph.renaming.GraphIdentifierLookupRenamer.ops._
import wdl.transforms.base.wdlom2wom.expression.renaming.IdentifierLookupRenamer.ops._

package object renaming {

  implicit val workflowDefinitionIdentifierRenamer: GraphIdentifierLookupRenamer[WorkflowDefinitionElement] =
    new GraphIdentifierLookupRenamer[WorkflowDefinitionElement] {
      override def renameIdentifiers(a: WorkflowDefinitionElement, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement],
        graphIdentifierLookupRenamer: GraphIdentifierLookupRenamer[WorkflowGraphElement]
      ): WorkflowDefinitionElement = {
        val newInputsSection = a.inputsSection map { inputSection =>
          inputSection.copy(inputDeclarations = inputSection.inputDeclarations.map { id =>
            val renamingMapWithoutThisValue = renamingMap.filterNot(_._1 == id.name)
            id.copy(expression = id.expression.map(_.renameIdentifiers(renamingMapWithoutThisValue)))
          })
        }

        val newGraphElements = a.graphElements map { ge => ge.renameIdentifiers(renamingMap) }

        val newOutputsSection = a.outputsSection map { outputSection =>
          outputSection.copy(outputs = outputSection.outputs.map { od =>
            od.copy(expression = od.expression.renameIdentifiers(renamingMap))
          })
        }

        a.copy(inputsSection = newInputsSection, graphElements = newGraphElements, outputsSection = newOutputsSection)
      }
    }

  implicit val graphElementIdentifierRenamer = new GraphIdentifierLookupRenamer[WorkflowGraphElement] {
    override def renameIdentifiers(a: WorkflowGraphElement, renamingMap: Map[String, String])(implicit
      expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement],
      graphIdentifierLookupRenamer: GraphIdentifierLookupRenamer[WorkflowGraphElement]
    ): WorkflowGraphElement = a match {
      case e: ScatterElement => e.renameIdentifiers(renamingMap)
      case e: IfElement => e.renameIdentifiers(renamingMap)
      case e: CallElement => e.renameIdentifiers(renamingMap)
      case e: InputDeclarationElement => e.renameIdentifiers(renamingMap)
      case e: IntermediateValueDeclarationElement => e.renameIdentifiers(renamingMap)
      case e: OutputDeclarationElement => e.renameIdentifiers(renamingMap)
    }
  }

  implicit val scatterElementIdentifierRenamer: GraphIdentifierLookupRenamer[ScatterElement] =
    new GraphIdentifierLookupRenamer[ScatterElement] {
      override def renameIdentifiers(a: ScatterElement, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement],
        graphIdentifierLookupRenamer: GraphIdentifierLookupRenamer[WorkflowGraphElement]
      ): ScatterElement =
        a.copy(
          scatterExpression = a.scatterExpression.renameIdentifiers(renamingMap),
          graphElements = a.graphElements.map(graphIdentifierLookupRenamer.renameIdentifiers(_, renamingMap))
        )
    }

  implicit val ifElementIdentifierRenamer: GraphIdentifierLookupRenamer[IfElement] =
    new GraphIdentifierLookupRenamer[IfElement] {
      override def renameIdentifiers(a: IfElement, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement],
        graphIdentifierLookupRenamer: GraphIdentifierLookupRenamer[WorkflowGraphElement]
      ): IfElement =
        a.copy(
          conditionExpression = a.conditionExpression.renameIdentifiers(renamingMap),
          graphElements = a.graphElements.map(graphIdentifierLookupRenamer.renameIdentifiers(_, renamingMap))
        )
    }

  implicit val callElementIdentifierRenamer: GraphIdentifierLookupRenamer[CallElement] =
    new GraphIdentifierLookupRenamer[CallElement] {
      override def renameIdentifiers(a: CallElement, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement],
        graphIdentifierLookupRenamer: GraphIdentifierLookupRenamer[WorkflowGraphElement]
      ): CallElement =
        a.copy(body = a.body map { bodyElement =>
          bodyElement.copy(inputs = bodyElement.inputs.map { callInput =>
            callInput.copy(value = callInput.value.renameIdentifiers(renamingMap))
          })
        })
    }

  implicit val inputDeclarationIdentifierRenamer: GraphIdentifierLookupRenamer[InputDeclarationElement] =
    new GraphIdentifierLookupRenamer[InputDeclarationElement] {
      override def renameIdentifiers(a: InputDeclarationElement, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement],
        graphIdentifierLookupRenamer: GraphIdentifierLookupRenamer[WorkflowGraphElement]
      ): InputDeclarationElement = {
        val renamingMapWithoutThisValue = renamingMap.filterNot(_._2 == a.name)
        a.copy(expression = a.expression.map(_.renameIdentifiers(renamingMapWithoutThisValue)))
      }
    }

  implicit val intermediateDeclarationIdentifierRenamer
    : GraphIdentifierLookupRenamer[IntermediateValueDeclarationElement] =
    new GraphIdentifierLookupRenamer[IntermediateValueDeclarationElement] {
      override def renameIdentifiers(a: IntermediateValueDeclarationElement, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement],
        graphIdentifierLookupRenamer: GraphIdentifierLookupRenamer[WorkflowGraphElement]
      ): IntermediateValueDeclarationElement = {
        val renamingMapWithoutThisValue = renamingMap.filterNot(_._2 == a.name)
        a.copy(expression = a.expression.renameIdentifiers(renamingMapWithoutThisValue))
      }
    }

  implicit val outputDeclarationIdentifierRenamer: GraphIdentifierLookupRenamer[OutputDeclarationElement] =
    new GraphIdentifierLookupRenamer[OutputDeclarationElement] {
      override def renameIdentifiers(a: OutputDeclarationElement, renamingMap: Map[String, String])(implicit
        expressionElementRenamer: IdentifierLookupRenamer[ExpressionElement],
        graphIdentifierLookupRenamer: GraphIdentifierLookupRenamer[WorkflowGraphElement]
      ): OutputDeclarationElement = {
        val renamingMapWithoutThisValue = renamingMap.filterNot(_._2 == a.name)
        a.copy(expression = a.expression.renameIdentifiers(renamingMapWithoutThisValue))
      }
    }
}
