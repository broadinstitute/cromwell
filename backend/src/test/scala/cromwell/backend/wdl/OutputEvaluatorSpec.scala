package cromwell.backend.wdl

import cats.data.{NonEmptyList, Validated}
import cats.syntax.validated._
import cromwell.backend.wdl.OutputEvaluator.{InvalidJobOutputs, JobOutputsEvaluationException, ValidJobOutputs}
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey}
import cromwell.core.NoIoFunctionSet
import cromwell.util.WomMocks
import lenthall.validation.ErrorOr.ErrorOr
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wom.JobOutput
import wom.callable.Callable.{InputDefinition, OutputDefinition, RequiredInputDefinition}
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.WomIdentifier
import wom.types.{WdlIntegerType, WdlType}
import wom.values.{WdlFile, WdlInteger, WdlValue}

class OutputEvaluatorSpec extends FlatSpec with Matchers with Mockito {
  behavior of "OutputEvaluator"
  
  // Depends on an input
  def o1Expression = new WomExpression {
    override def sourceString: String = "o1"
    override def inputs: Set[String] = Set("input")
    override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
      Validated.fromOption(inputValues.get("input"), NonEmptyList.one("Can't find a value for 'input'"))
    }
    override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = ???
    override def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] = ???
  }

  // Depends on a previous output
  def o2Expression = new WomExpression {
    override def sourceString: String = "o2"
    override def inputs: Set[String] = Set("o1")
    override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
      Validated.fromOption(inputValues.get("o1"), NonEmptyList.one("Can't find a value for 'o1'"))
    }
    override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = ???
    override def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] = ???
  }

  def invalidWomExpression1 = new WomExpression {
    override def sourceString: String = "invalid1"
    override def inputs: Set[String] = Set.empty
    override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
      "Invalid expression 1".invalidNel
    }
    override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = {
      "Invalid expression 1".invalidNel
    }
    override def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] = {
      "Invalid expression 1".invalidNel
    }
  }

  def invalidWomExpression2 = new WomExpression {
    override def sourceString: String = "invalid2"
    override def inputs: Set[String] = Set.empty
    override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
      "Invalid expression 2".invalidNel
    }
    override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = {
      "Invalid expression 2".invalidNel
    }
    override def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] = {
      "Invalid expression 2".invalidNel
    }
  }
  
  val exception = new Exception("Expression evaluation exception")
  
  def throwingWomExpression = new WomExpression {
    override def sourceString: String = "throwing"
    override def inputs: Set[String] = Set.empty
    override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
      throw exception
    }
    override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = {
      throw exception
    }
    override def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] = {
      throw exception
    }
  }

  val mockInputs: Map[InputDefinition, WdlValue] = Map(
    RequiredInputDefinition("input", WdlIntegerType) -> WdlInteger(5)
  )
  
  it should "evaluate valid jobs outputs" in {
    val mockOutputs = List (
      OutputDefinition("o1", WdlIntegerType, o1Expression),
      OutputDefinition("o2", WdlIntegerType, o2Expression)
    )
    
    val call = WomMocks.mockTaskCall(WomIdentifier("call"), WomMocks.EmptyTaskDefinition.copy(outputs = mockOutputs))
    val key = BackendJobDescriptorKey(call, None, 1)
    val jobDescriptor = BackendJobDescriptor(null, key, null, mockInputs, null, null)
    
    OutputEvaluator.evaluateOutputs(jobDescriptor, NoIoFunctionSet) match {
      case ValidJobOutputs(outputs) => outputs shouldBe Map(
        "o1" -> JobOutput(WdlInteger(5)),
        "o2" -> JobOutput(WdlInteger(5))
      )
      case _ => fail("Failed to evaluate outputs")
    }
  }

  it should "return an InvalidJobOutputs if the evaluation returns ErrorOrs" in {
    val mockOutputs = List (
      OutputDefinition("o1", WdlIntegerType, o1Expression),
      OutputDefinition("invalid1", WdlIntegerType, invalidWomExpression1),
      OutputDefinition("invalid2", WdlIntegerType, invalidWomExpression2)
    )

    val call = WomMocks.mockTaskCall(WomIdentifier("call"), WomMocks.EmptyTaskDefinition.copy(outputs = mockOutputs))
    val key = BackendJobDescriptorKey(call, None, 1)
    val jobDescriptor = BackendJobDescriptor(null, key, null, mockInputs, null, null)

    OutputEvaluator.evaluateOutputs(jobDescriptor, NoIoFunctionSet) match {
      case InvalidJobOutputs(errors) => errors shouldBe NonEmptyList.of(
        "Invalid expression 1", "Invalid expression 2"
      )
      case _ => fail("Output evaluation should have failed")
    }
  }

  it should "return an JobOutputsEvaluationException if the evaluation throws an exception" in {
    val mockOutputs = List (
      OutputDefinition("o1", WdlIntegerType, o1Expression),
      OutputDefinition("invalid1", WdlIntegerType, throwingWomExpression)
    )

    val call = WomMocks.mockTaskCall(WomIdentifier("call"), WomMocks.EmptyTaskDefinition.copy(outputs = mockOutputs))
    val key = BackendJobDescriptorKey(call, None, 1)
    val jobDescriptor = BackendJobDescriptor(null, key, null, mockInputs, null, null)

    OutputEvaluator.evaluateOutputs(jobDescriptor, NoIoFunctionSet) match {
      case JobOutputsEvaluationException(e) => e shouldBe exception
      case _ => fail("Output evaluation should have failed")
    }
  }
}
