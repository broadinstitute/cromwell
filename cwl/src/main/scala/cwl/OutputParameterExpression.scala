package cwl

import common.validation.ErrorOr.ErrorOr
import common.validation.IOChecked.{IOChecked, _}
import cwl.CwlType.CwlType
import shapeless.Poly1
import wom.expression.{EmptyIoFunctionSet, FileEvaluation, IoFunctionSet}
import wom.types._
import wom.values.WomValue

import scala.Function.const
import scala.concurrent.{ExecutionContext, Future}

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
                                    cwlExpressionType: WomType): IOChecked[WomValue] = {
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
                                   )(outputBinding: CommandOutputBinding): IOChecked[Set[FileEvaluation]] = {
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

    fromOutputBinding.orElse(fromType).getOrElse(s"Cannot evaluate ${parameter.toString}".invalidIOChecked).toErrorOr
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
  override def evaluateFiles(inputs: Map[String, WomValue], unused: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[FileEvaluation]] = {
    import cats.syntax.apply._

    // Ignore the supplied ioFunctionSet and use a custom stubbed IoFunctionSet. This is better than a real I/O function set
    // because in the context of file evaluation we don't care about the results of these operations. The NoIoFunctionSet
    // that otherwise would be used throws for all of its operations which doesn't fly for the way our CWL evaluation works.
    val stubbedIoFunctionSet = new EmptyIoFunctionSet {
      override def size(path: String): Future[Long] = Future.successful(0L)
      override def isDirectory(path: String): Future[Boolean] = Future.successful(false)
      override def ec: ExecutionContext = scala.concurrent.ExecutionContext.global
    }
    def fromOutputBinding: IOChecked[Set[FileEvaluation]] = parameter
      .outputBinding
      .map(evaluateOutputBindingFiles(inputs, stubbedIoFunctionSet, parameter.secondaryFiles, coerceTo))
      .getOrElse(Set.empty[FileEvaluation].validIOChecked)

    def fromType: IOChecked[Set[FileEvaluation]] = parameter
      .`type`
      .map(_.fold(MyriadOutputTypeToWomFiles).apply(evaluateOutputBindingFiles(inputs, stubbedIoFunctionSet, parameter.secondaryFiles, coerceTo)))
      .getOrElse(Set.empty[FileEvaluation].validIOChecked)

    val optional: Boolean = parameter.`type`.exists(_.fold(OutputTypeIsOptional))

    ((fromOutputBinding, fromType) mapN (_ ++ _) map { _ map { _.copy(optional = optional) } }).toErrorOr
  }
}

object OutputTypeIsOptional extends Poly1 {
  implicit val one: Case.Aux[MyriadOutputInnerType, Boolean] = at[MyriadOutputInnerType] { const(false) }

  implicit val arr: Case.Aux[Array[MyriadOutputInnerType], Boolean] = at[Array[MyriadOutputInnerType]] {
    // Possibly too broad, would return true for just single 'null'.
    _.exists(_.select[CwlType].contains(CwlType.Null))
  }
}
