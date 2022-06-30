package cromwell.backend.impl.bcs

import common.validation.ErrorOr.ErrorOr
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.standard.{StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.filesystems.oss.OssPathBuilder
import cromwell.filesystems.oss.OssPathBuilder.{InvalidOssPath, PossiblyValidRelativeOssPath, ValidFullOssPath}
import wom.graph.CommandCallNode
import wom.values.WomGlobFile

import scala.concurrent.Future

final case class BcsExpressionFunctions(override val standardParams: StandardExpressionFunctionsParams)
  extends StandardExpressionFunctions(standardParams) {

  override def preMapping(str: String) = {
    OssPathBuilder.validateOssPath(str) match {
      case _: ValidFullOssPath => str
      case PossiblyValidRelativeOssPath => callContext.root.resolve(str.stripPrefix("/")).pathAsString
      case invalid: InvalidOssPath => throw new IllegalArgumentException(invalid.errorMessage)
    }
  }

  // TODO: BCS: When globs are supported this override should be removed.
  // https://github.com/broadinstitute/cromwell/issues/3519
  // See usages of cwl.CommandLineTool.CwlOutputJson
  override def glob(pattern: String): Future[Seq[String]] = {
    if (pattern == "cwl.output.json") {
      Future.successful(Nil)
    } else {
      super.glob(pattern)
    }
  }

  // TODO: BCS: When globs are supported this override should be removed.
  // https://github.com/broadinstitute/cromwell/issues/3519
  // See usages of cwl.CommandLineTool.CwlOutputJson
  override def findGlobOutputs(call: CommandCallNode, jobDescriptor: BackendJobDescriptor): ErrorOr[List[WomGlobFile]] = {
    val base = super.findGlobOutputs(call, jobDescriptor)
    base.map(_.filterNot(_.value == "cwl.output.json"))
  }
}
