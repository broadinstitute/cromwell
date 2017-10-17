package wom.expression

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import lenthall.validation.ErrorOr.ErrorOr
import wdl.types.WdlType
import wdl.values.{WdlFile, WdlFloat, WdlValue}
import wom.graph.{GraphNode, InstantiatedExpression}
import wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort}
import cats.syntax.validated._
import cats.instances.list._
import cats.syntax.traverse._

import scala.concurrent.Future
import scala.util.Try

trait WomExpression {
  def sourceString: String
  def inputs: Set[String]
  def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue]
  def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType]
  def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]]
  def linkWithInputs(graphNodeSetter: GraphNode.GraphNodeSetter, inputMapping: Map[String, OutputPort]): ErrorOr[InstantiatedExpression] = {

    def linkInput(input: String): ErrorOr[(String, InputPort)] = if (inputMapping.contains(input)) {
      val upstreamPort = inputMapping(input)
      Valid((input, ConnectedInputPort(input, upstreamPort.womType, upstreamPort, graphNodeSetter.get)))
    } else {
      s"Expression cannot be connected without the input $input (provided: ${inputMapping.toString})".invalidNel
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

final case class PlaceholderWomExpression(inputs: Set[String], fixedWomType: WdlType) extends WomExpression {
  override def sourceString: String = "placeholder"
  override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] =
    Invalid(NonEmptyList.one(s"couldn't evaluate value from inputs $inputs\tfixedWomType\t$fixedWomType\tinputValues\t$inputValues"))
  override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] =
    Valid(fixedWomType)
  override def evaluateFiles(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] =
    Valid(Set.empty)
}

// TODO: Flesh this out (https://github.com/broadinstitute/cromwell/issues/2521)
trait IoFunctionSet {
  def readFile(path: String): Future[String]
  def writeFile(path: String, content: String): Future[WdlFile]
  def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile]
  def glob(path: String, pattern: String): Seq[String]
  def size(params: Seq[Try[WdlValue]]): Try[WdlFloat]
}

case object PlaceholderIoFunctionSet extends IoFunctionSet {
  override def readFile(path: String): Future[String] = ???
  override def writeFile(path: String, content: String): Future[WdlFile] = ???
  override def stdout(params: Seq[Try[WdlValue]]) = ???
  override def stderr(params: Seq[Try[WdlValue]]) = ???
  override def glob(path: String, pattern: String) = ???
  override def size(params: Seq[Try[WdlValue]]) = ???
}
