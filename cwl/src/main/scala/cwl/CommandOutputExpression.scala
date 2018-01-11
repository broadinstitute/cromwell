package cwl

import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.types._
import wom.values.{WomFile, WomValue}

case class CommandOutputExpression(outputBinding: CommandOutputBinding,
                                   override val cwlExpressionType: WomType,
                                   override val inputs: Set[String],
                                   secondaryFilesOption: Option[SecondaryFiles] = None,
                                   formatOption: Option[StringOrExpression] = None //only valid when type: File
                                  ) extends CwlWomExpression {

  // TODO WOM: outputBinding.toString is probably not the best representation of the expression source
  override def sourceString = outputBinding.toString

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
    CommandOutputBinding.generateOutputWomValue(
      inputValues,
      ioFunctionSet,
      cwlExpressionType,
      outputBinding,
      secondaryFilesOption,
      formatOption
    )
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
  override def evaluateFiles(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] ={
    CommandOutputBinding.getOutputWomFiles(
      inputs,
      coerceTo,
      outputBinding,
      secondaryFilesOption
    )
  }
}
