package cwl

import cats.syntax.option._
import cats.syntax.validated._
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._

final case class WorkflowStepInputExpression(input: WorkflowStepInput, override val cwlExpressionType: WomType, graphInputs: Set[String], override val expressionLib: ExpressionLib)(implicit parentName: ParentName) extends CwlWomExpression {

  override def sourceString = input.toString

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) = {
    (input.valueFrom, input.source) match {
      case (None, Some(WorkflowStepInputSource.String(id))) =>
        inputValues.
          get(FullyQualifiedName(id).id).
          toValidNel(s"could not find id $id in typeMap\n${inputValues.mkString("\n")}\nwhen evaluating $input.  Graph Inputs were ${graphInputs.mkString("\n")}")
      case _ => s"Could not do evaluateValue(${input.valueFrom}, ${input.source}), most likely it has not been implemented yet".invalidNel
    }
  }

  //this is input, so not producing any output files
  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) =
    Set.empty[WomFile].validNel

  override def inputs = graphInputs ++ input.source.toSet.flatMap{ inputSource: InputSource => inputSource match {
    case WorkflowStepInputSource.String(s) => Set(FullyQualifiedName(s).id)
    case WorkflowStepInputSource.StringArray(sa) => sa.map(FullyQualifiedName(_).id).toSet
  }}
}
