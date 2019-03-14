package cromwell.backend

import java.util.concurrent.Executors

import cats.data.{NonEmptyList, Validated}
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.OutputEvaluator.{InvalidJobOutputs, JobOutputsEvaluationException, ValidJobOutputs}
import cromwell.core.CallOutputs
import cromwell.util.WomMocks
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import wom.callable.Callable.{InputDefinition, OutputDefinition, RequiredInputDefinition}
import wom.expression.{FileEvaluation, IoFunctionSet, NoIoFunctionSet, WomExpression}
import wom.graph.WomIdentifier
import wom.types.{WomIntegerType, WomType}
import wom.values.{WomInteger, WomValue}

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, ExecutionContextExecutor}

class OutputEvaluatorSpec extends FlatSpec with Matchers with Mockito {
  behavior of "OutputEvaluator"

  val FutureTimeout = 20.seconds
  final implicit val blockingEc: ExecutionContextExecutor = ExecutionContext.fromExecutor(
    Executors.newCachedThreadPool()
  )
  
  // Depends on an input
  def o1Expression = new WomExpression {
    override def sourceString: String = "o1"
    override def inputs: Set[String] = Set("input")
    override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      Validated.fromOption(inputValues.get("input"), NonEmptyList.one("Can't find a value for 'input'"))
    }
    override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = throw new UnsupportedOperationException
    override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[FileEvaluation]] = throw new UnsupportedOperationException
  }

  // Depends on a previous output
  def o2Expression = new WomExpression {
    override def sourceString: String = "o2"
    override def inputs: Set[String] = Set("o1")
    override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      Validated.fromOption(inputValues.get("o1"), NonEmptyList.one("Can't find a value for 'o1'"))
    }
    override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = throw new UnsupportedOperationException
    override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[FileEvaluation]] = throw new UnsupportedOperationException
  }

  def invalidWomExpression1 = new WomExpression {
    override def sourceString: String = "invalid1"
    override def inputs: Set[String] = Set.empty
    override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      "Invalid expression 1".invalidNel
    }
    override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = {
      "Invalid expression 1".invalidNel
    }
    override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[FileEvaluation]] = {
      "Invalid expression 1".invalidNel
    }
  }

  def invalidWomExpression2 = new WomExpression {
    override def sourceString: String = "invalid2"
    override def inputs: Set[String] = Set.empty
    override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      "Invalid expression 2".invalidNel
    }
    override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = {
      "Invalid expression 2".invalidNel
    }
    override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[FileEvaluation]] = {
      "Invalid expression 2".invalidNel
    }
  }
  
  val exception = new Exception("Expression evaluation exception")
  
  def throwingWomExpression = new WomExpression {
    override def sourceString: String = "throwing"
    override def inputs: Set[String] = Set.empty
    override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      throw exception
    }
    override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = {
      throw exception
    }
    override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[FileEvaluation]] = {
      throw exception
    }
  }

  val mockInputs: Map[InputDefinition, WomValue] = Map(
    RequiredInputDefinition("input", WomIntegerType) -> WomInteger(5)
  )
  
  it should "evaluate valid jobs outputs" in {
    val mockOutputs = List (
      OutputDefinition("o1", WomIntegerType, o1Expression),
      OutputDefinition("o2", WomIntegerType, o2Expression)
    )
    
    val call = WomMocks.mockTaskCall(WomIdentifier("call"), WomMocks.EmptyTaskDefinition.copy(outputs = mockOutputs))
    val key = BackendJobDescriptorKey(call, None, 1)
    val jobDescriptor = BackendJobDescriptor(null, key, null, mockInputs, null, None, null)
    
    Await.result(OutputEvaluator.evaluateOutputs(jobDescriptor, NoIoFunctionSet), FutureTimeout) match {
      case ValidJobOutputs(outputs) => outputs shouldBe CallOutputs(Map(
        jobDescriptor.taskCall.outputPorts.find(_.name == "o1").get -> WomInteger(5),
        jobDescriptor.taskCall.outputPorts.find(_.name == "o2").get -> WomInteger(5)
      ))
      case _ => fail("Failed to evaluate outputs")
    }
  }

  it should "return an InvalidJobOutputs if the evaluation returns ErrorOrs" in {
    val mockOutputs = List (
      OutputDefinition("o1", WomIntegerType, o1Expression),
      OutputDefinition("invalid1", WomIntegerType, invalidWomExpression1),
      OutputDefinition("invalid2", WomIntegerType, invalidWomExpression2)
    )

    val call = WomMocks.mockTaskCall(WomIdentifier("call"), WomMocks.EmptyTaskDefinition.copy(outputs = mockOutputs))
    val key = BackendJobDescriptorKey(call, None, 1)
    val jobDescriptor = BackendJobDescriptor(null, key, null, mockInputs, null, None, null)

    Await.result(OutputEvaluator.evaluateOutputs(jobDescriptor, NoIoFunctionSet), FutureTimeout) match {
      case InvalidJobOutputs(errors) => errors shouldBe NonEmptyList.of(
        "Bad output 'invalid1': Invalid expression 1", "Bad output 'invalid2': Invalid expression 2"
      )
      case _ => fail("Output evaluation should have failed")
    }
  }

  it should "return an JobOutputsEvaluationException if the evaluation throws an exception" in {
    val mockOutputs = List (
      OutputDefinition("o1", WomIntegerType, o1Expression),
      OutputDefinition("invalid1", WomIntegerType, throwingWomExpression)
    )

    val call = WomMocks.mockTaskCall(WomIdentifier("call"), WomMocks.EmptyTaskDefinition.copy(outputs = mockOutputs))
    val key = BackendJobDescriptorKey(call, None, 1)
    val jobDescriptor = BackendJobDescriptor(null, key, null, mockInputs, null, None, null)

    Await.result(OutputEvaluator.evaluateOutputs(jobDescriptor, NoIoFunctionSet), FutureTimeout) match {
      case JobOutputsEvaluationException(e) => e shouldBe exception
      case _ => fail("Output evaluation should have failed")
    }
  }
}
