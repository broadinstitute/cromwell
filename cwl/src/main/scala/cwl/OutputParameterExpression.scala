package cwl

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cwl.CwlType.CwlType
import shapeless.Poly1
import wom.expression.{FileEvaluation, IoFunctionSet}
import wom.types._
import wom.values.{WomFile, WomValue}

import scala.Function.const

case class OutputParameterExpression(parameter: OutputParameter,
                                     override val cwlExpressionType: WomType,
                                     override val inputs: Set[String],
                                     override val expressionLib: ExpressionLib,
                                     schemaDefRequirement: SchemaDefRequirement) extends CwlWomExpression {

  override def sourceString = parameter.toString

  override def cacheString: String = parameter.cacheString

  private def evaluateOutputBinding(inputValues: Map[String, WomValue],
                                    ioFunctionSet: IoFunctionSet,
                                    secondaryFilesOption: Option[SecondaryFiles],
                                    formatOption: Option[StringOrExpression]
                                   )(outputBinding: CommandOutputBinding,
                                    cwlExpressionType: WomType) = {
    CommandOutputBinding.generateOutputWomValue(
      inputValues,
      ioFunctionSet,
      cwlExpressionType,
      outputBinding,
      secondaryFilesOption,
      formatOption,
      expressionLib
    )
  }

  private def evaluateOutputBindingFiles(inputValues: Map[String, WomValue],
                                    ioFunctionSet: IoFunctionSet,
                                    secondaryFilesOption: Option[SecondaryFiles],
                                    coerceTo: WomType
                                   )(outputBinding: CommandOutputBinding) = {
    CommandOutputBinding.getOutputWomFiles(
      inputValues,
      coerceTo,
      outputBinding,
      secondaryFilesOption,
      ioFunctionSet,
      expressionLib
    )
  }
  
  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
    def fromOutputBinding =
      parameter.outputBinding.map(evaluateOutputBinding(inputValues, ioFunctionSet, parameter.secondaryFiles, parameter.format)(_, cwlExpressionType))

    def fromType =
      parameter.`type`.map(_.fold(MyriadOutputTypeToWomValue).apply(
        evaluateOutputBinding(inputValues, ioFunctionSet, parameter.secondaryFiles, parameter.format),
        schemaDefRequirement
      ))

    fromOutputBinding.orElse(fromType).getOrElse(s"Cannot evaluate ${parameter.toString}".invalidNel)
  }

  /**
    * Returns the list of files that _will be_ output after the command is run, unless they are optional and if so _may be_.
    *
    * In CWL, a list of outputs is specified as glob, say `*.bam`, plus a list of secondary files that may be in the
    * form of paths or specified using carets such as `^.bai`.
    *
    * The coerceTo may be one of four different values:
    * - WomMaybePopulatedFileType
    * - WomArrayType(WomMaybePopulatedFileType)
    * - WomMaybeListedDirectoryType
    * - WomArrayType(WomMaybeListedDirectoryType) (Possible according to the way the spec is written, but not likely?)
    */
  override def evaluateFiles(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[FileEvaluation]] = {
    import cats.syntax.apply._
    
    def fromOutputBinding: ErrorOr[Set[WomFile]] = parameter
      .outputBinding
      .map(evaluateOutputBindingFiles(inputs, ioFunctionSet, parameter.secondaryFiles, coerceTo))
      .getOrElse(Set.empty[WomFile].validNel)

    def fromType: ErrorOr[Set[WomFile]] = parameter
      .`type`
      .map(_.fold(MyriadOutputTypeToWomFiles).apply(evaluateOutputBindingFiles(inputs, ioFunctionSet, parameter.secondaryFiles, coerceTo)))
      .getOrElse(Set.empty[WomFile].validNel)

    val optional: Boolean = parameter.`type`.exists(_.fold(OutputTypeIsOptional))

    (fromOutputBinding, fromType) mapN (_ ++ _) map { _ map { FileEvaluation(_, optional) } }
  }
}

object OutputTypeIsOptional extends Poly1 {
  implicit val one: Case.Aux[MyriadOutputInnerType, Boolean] = at[MyriadOutputInnerType] { const(false) }

  implicit val arr: Case.Aux[Array[MyriadOutputInnerType], Boolean] = at[Array[MyriadOutputInnerType]] {
    // Possibly too broad, would return true for just single 'null'.
    _.exists(_.select[CwlType].contains(CwlType.Null))
  }
}
