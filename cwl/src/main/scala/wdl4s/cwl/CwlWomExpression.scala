package wdl4s.cwl

import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wdl.types.WdlType
import wdl4s.wdl.values.{WdlFile, WdlGlobFile, WdlInteger, WdlSingleFile, WdlValue}
import wdl4s.wom.expression.{IoFunctionSet, WomExpression}

import scala.concurrent.Await
import scala.concurrent.duration.Duration

object CwlWomExpression {
  // FIXME This implementation is obviously a joke only designed as a placeholder for real expression evaluation
  // so that we get 3step.wdl and 3step.cwl running
  def apply(commandOutputParameter: CommandOutputParameter, wdlType: WdlType): WomExpression = {
    val globFileName = for {
      outputBinding <- commandOutputParameter.outputBinding
      glob <- outputBinding.glob
      fileName <- glob.select[String]
    } yield fileName

    val outputEval = for {
      outputBinding <- commandOutputParameter.outputBinding
      outputEval <- outputBinding.outputEval
      expr <- outputEval.select[String]
    } yield expr

    new WomExpression {
      def evaluateGlob(globName: String, ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
        WdlSingleFile(ioFunctionSet.glob("", globName).head).validNel
      }

      def evaluateOutputEval(expr: String, glob: WdlValue, ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
        if (expr == "$(self[0].contents.toInt)") {
          val content = Await.result(ioFunctionSet.readFile(glob.valueString), Duration.Inf)
          WdlInteger(content.trim.toInt).validNel
        } else s"Can't evaluate $expr".invalidNel
      }

      override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
        val globResult = globFileName map {fileName => evaluateGlob(fileName, ioFunctionSet) }
        outputEval.map(evaluateOutputEval(_, globResult.get.getOrElse(???), ioFunctionSet)).getOrElse(globResult.get)
      }

      override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = ???

      override def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] = {
        globFileName.map(WdlGlobFile).toSet[WdlFile].validNel
      }
      override def inputs: Set[String] = ???
    }
  }

  def apply(expression: String): WomExpression = {
    new WomExpression {
      override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
        expression.split('.').toList match {
          case "inputs" :: value :: Nil => inputValues(value).validNel
          case _ => s"Can't find value for $expression".invalidNel
        }

      }
      override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = ???
      override def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] = ???
      override def inputs: Set[String] = ???
    }
  }
}
