package cwl

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import cats.syntax.option._
import cwl.WorkflowStepInput.InputSource
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.Validation._
import shapeless.{Inl, Poly1}
import wdl.types._
import wdl.values.{WdlArray, WdlFile, WdlGlobFile, WdlMap, WdlString, WdlValue}
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort}
import wom.graph.{GraphNode, InstantiatedExpression}
import cats.syntax.validated._
import cats.instances.list._
import cats.syntax.traverse._

sealed trait CwlWomExpression extends WomExpression {

  def cwlExpressionType: WdlType

  override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = cwlExpressionType.validNel

  override def linkWithInputs(graphNodeSetter: GraphNode.GraphNodeSetter, fullyQualifiedInputMapping: Map[String, OutputPort]): ErrorOr[InstantiatedExpression] = {

//    val inputMapping = fullyQualifiedInputMapping.map { case (key, value) => FullyQualifiedName(key).id -> value}
    val inputMapping = fullyQualifiedInputMapping

    println(s"input mpaping is $inputMapping")


    def linkInput(input: String): ErrorOr[(String, InputPort)] = if (inputMapping.contains(input)) {
      val upstreamPort = inputMapping(input)
      Valid((input, ConnectedInputPort(input, upstreamPort.womType, upstreamPort, graphNodeSetter.get)))
    } else {
      s"Expression cannot be connected without the input $input (provided:\n${inputMapping.mkString("\n")})".invalidNel
    }

    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap
    for {
      linkedInputList <- inputs.toList traverse linkInput
      linkedInputs = linkedInputList.toMap
      inputTypes = linkedInputs map { case (k, v) => k -> v.womType }
      evaluatedType <- evaluateType(inputTypes)
    } yield new InstantiatedExpression(this, evaluatedType, linkedInputs)
  }
}

case class CommandOutputExpression(outputBinding: CommandOutputBinding,
                                   override val cwlExpressionType: WdlType, override val inputs: Set[String]) extends CwlWomExpression {

  // TODO WOM: outputBinding.toString is probably not be the best representation of the outputBinding
  override def sourceString = outputBinding.toString
  override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
    val parameterContext = ParameterContext.Empty.withInputs(inputValues, ioFunctionSet)

    val wdlValue: WdlValue = outputBinding.commandOutputBindingToWdlValue(parameterContext, ioFunctionSet)
    val extractFile: WdlValue =
      wdlValue match {
        case WdlArray(WdlMaybeEmptyArrayType(WdlMapType(WdlStringType, WdlStringType)), Seq(WdlMap(WdlMapType(WdlStringType, WdlStringType), map))) =>
          map(WdlString("location"))
        case other => other
      }
    cwlExpressionType.coerceRawValue(extractFile).toErrorOr
  }

  /*
  TODO:
   DB: It doesn't make sense to me that this function returns type WdlFile but accepts a type to which it coerces.
   Wouldn't coerceTo always == WdlFileType, and if not then what?
   */
  override def evaluateFiles(inputs: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] ={

    val pc = ParameterContext.Empty.withInputs(inputs, ioFunctionSet)

    outputBinding.glob.toList.flatMap { globValue =>
      GlobEvaluator.globPaths(globValue, pc, ioFunctionSet).toList
    }.map{s:String => WdlGlobFile(s): WdlFile}.toSet.validNel[String]
  }
}

case class WorkflowStepInputExpression(input: WorkflowStepInput, override val cwlExpressionType: WdlType, val graphInputs: Set[String]) extends CwlWomExpression {

  override def sourceString = input.toString

  override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet) = {
    (input.valueFrom, input.source) match {
      case (None, Some(Inl(id: String))) =>
        inputValues.
          get(id).
          toValidNel(s"could not find id $id in typeMap ${inputValues.foreach(println)}\twhen evaluating $input")
      case _ => Invalid(NonEmptyList.one("could not decipher evaluateValue, most likely has not been implemented yet"))
    }
  }

  override def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType) = ???

  object InputSourceToFileNames extends Poly1{

    implicit def string = at[String]{s => Set(FullyQualifiedName(s).id)}

    implicit def array = at[Array[String]]{_.map(FullyQualifiedName(_).id).toSet}
  }

  override def inputs = graphInputs ++ input.source.toSet.flatMap{(_:InputSource).fold(InputSourceToFileNames)}
}


