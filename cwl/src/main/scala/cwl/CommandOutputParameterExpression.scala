package cwl

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.types._
import wom.values.{WomFile, WomValue}

case class CommandOutputParameterExpression(parameter: OutputParameter,
                                            override val cwlExpressionType: WomType,
                                            override val inputs: Set[String],
                                            override val expressionLib: ExpressionLib) extends CwlWomExpression {

  override def sourceString = parameter.toString

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
      expressionLib
    )
  }
  
  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
    def fromOutputBinding =
      parameter.outputBinding.map(evaluateOutputBinding(inputValues, ioFunctionSet, parameter.secondaryFiles, parameter.format)(_, cwlExpressionType))

    def fromType =
      parameter.`type`.map(_.fold(MyriadOutputTypeToWomValue).apply(evaluateOutputBinding(inputValues, ioFunctionSet, parameter.secondaryFiles, parameter.format)))

    fromOutputBinding.orElse(fromType).getOrElse(s"Cannot evaluate ${parameter.toString}".invalidNel)
  }

  /**
    * Returns the list of files that _will be_ output after the command is run.
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
  override def evaluateFiles(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] = {
    import cats.syntax.apply._
    
    def fromOutputBinding: ErrorOr[Set[WomFile]] = parameter
      .outputBinding
      .map(evaluateOutputBindingFiles(inputs, ioFunctionSet, parameter.secondaryFiles, coerceTo))
      .getOrElse(Set.empty[WomFile].validNel)

    def fromType: ErrorOr[Set[WomFile]] = parameter
      .`type`
      .map(_.fold(MyriadOutputTypeToWomFiles).apply(evaluateOutputBindingFiles(inputs, ioFunctionSet, parameter.secondaryFiles, coerceTo)))
      .getOrElse(Set.empty[WomFile].validNel)
    
    (fromOutputBinding, fromType) mapN (_ ++ _)
  }
}
