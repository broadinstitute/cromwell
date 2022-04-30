import cats.data.Validated._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers._
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.elements.{InputDeclarationElement, InputsSectionElement, IntermediateValueDeclarationElement, OutputDeclarationElement, OutputsSectionElement, PrimitiveTypeElement}
import wom.types.{WomSingleFileType, WomStringType}
import wdl.transforms.base.ast2wdlom.AstToWorkflowDefinitionElement._


class AstToWorkflowDefinitionElementSpec extends AnyFlatSpec{
  behavior of "Check Stdouts and Stderrs"

  val mockInputSectionStdout = Option(InputsSectionElement(Vector(InputDeclarationElement(PrimitiveTypeElement(WomSingleFileType), "i", Some(StdoutElement)))))
  val mockInputSectionStderr = Option(InputsSectionElement(Vector(InputDeclarationElement(PrimitiveTypeElement(WomSingleFileType), "i", Some(StderrElement)))))
  val mockInputSectionNonStd = Option(InputsSectionElement(Vector(InputDeclarationElement(PrimitiveTypeElement(WomStringType), "more", Some(StringLiteral("more"))))))


  val mockIntermediatesStdout = Vector(IntermediateValueDeclarationElement(PrimitiveTypeElement(WomSingleFileType), "y", StdoutElement))
  val mockIntermediatesStderr = Vector(IntermediateValueDeclarationElement(PrimitiveTypeElement(WomSingleFileType), "y", StderrElement))
  val mockIntermediatesNonStd = Vector(IntermediateValueDeclarationElement(PrimitiveTypeElement(WomStringType), "here", StringLiteral("here")))

  val mockOutputSectionStdout = Option(OutputsSectionElement(Vector(OutputDeclarationElement(PrimitiveTypeElement(WomSingleFileType), "s", StdoutElement))))
  val mockOutputSectionStderr = Option(OutputsSectionElement(Vector(OutputDeclarationElement(PrimitiveTypeElement(WomSingleFileType), "s", StderrElement))))
  val mockOutputSectionNonStd = Option(OutputsSectionElement(Vector(OutputDeclarationElement(PrimitiveTypeElement(WomStringType), "more", StringLiteral("more")))))


  it should "return an error when there is an stdout in input section" in {
    val testInputs = checkDisallowedInputElement(mockInputSectionStdout, StdoutElement, "stdout")
    testInputs match {
      case Valid(_) => fail("Input section contained stdout. Should have failed.")
      case Invalid(e) => e.head should be("Workflow cannot have stdout expression in input section at workflow-level.")
    }
  }

  it should "return an error when there is an stderr in input section" in {
    val testInputs = checkDisallowedInputElement(mockInputSectionStderr, StderrElement, "stderr")
    testInputs match {
      case Valid(_) => fail("Input section contained stderr. Should have failed.")
      case Invalid(e) => e.head should be("Workflow cannot have stderr expression in input section at workflow-level.")
    }
  }

  it should "not return an error for non-stdout/stderr in the inputs section" in {
    val testInputs = checkDisallowedInputElement(mockInputSectionNonStd, StdoutElement, "non-stdout/stderr")
    testInputs match {
      case Valid(_) => // No action
      case Invalid(_) => fail("Check shouldn't have returned error as input section had non-stdout/stderr inputs")
    }
  }

  it should "return an error when there is an stdout in output section" in {
    val testOutputs = checkDisallowedOutputElement(mockOutputSectionStdout, StdoutElement, "stdout")
    testOutputs match {
      case Valid(_) => fail("Output section contained stdout. Should have failed.")
      case Invalid(e) => e.head should be("Workflow cannot have stdout expression in output section at workflow-level.")
    }
  }

  it should "return an error when there is an stderr in output section" in {
    val testOutputs = checkDisallowedOutputElement(mockOutputSectionStderr, StderrElement, "stderr")
    testOutputs match {
      case Valid(_) => fail("Output section contained stderr. Should have failed.")
      case Invalid(e) => e.head should be("Workflow cannot have stderr expression in output section at workflow-level.")
    }
  }
  it should "not return an error for non-stdout/stderr in the outputs section" in {
    val testOutputs = checkDisallowedOutputElement(mockOutputSectionNonStd, StdoutElement, "non-stdout/stderr")
    testOutputs match {
      case Valid(_) => // No action
      case Invalid(_) => fail("Check shouldn't have returned error as output section had non-stdout/stderr outputs")
    }
  }

  it should "return an error when there is an stdout at intermediate declaration section" in {
    val testIntermediates = checkDisallowedIntermediates(mockIntermediatesStdout, StdoutElement, "stdout")
    testIntermediates match {
      case Valid(_) => fail("Intermediate section contained stdout. Should have failed.")
      case Invalid(e) => e.head should be("Workflow cannot have stdout expression at intermediate declaration section at workflow-level.")
    }
  }
  it should "return an error when there is an stderr at intermediate declaration section" in {
    val testIntermediates = checkDisallowedIntermediates(mockIntermediatesStderr, StderrElement, "stderr")
    testIntermediates match {
      case Valid(_) => fail("Intermediate section contained stderr. Should have failed.")
      case Invalid(e) => e.head should be("Workflow cannot have stderr expression at intermediate declaration section at workflow-level.")
    }
  }

  it should "not return an error for non-stdout/stderr in the intermediates section" in {
    val testIntermediates = checkDisallowedIntermediates(mockIntermediatesNonStd, StdoutElement, "non-stdout/stderr")
    testIntermediates match {
      case Valid(_) => // No action
      case Invalid(_) => fail("Check shouldn't have returned error as intermediate section had non-stdout/stderr intermediates.")
    }
  }

}