package cwl

import shapeless.{Path => _, _}
import cats.syntax.traverse._
import cats.instances.list._
import cats.effect.IO
import common.validation.Parse._
import cats.{Applicative, Monad}
import cats.data.EitherT._
import cats.data.EitherT
import common.validation.ErrorOr._

object AddEmbeddedCwl extends Poly1 {

  implicit val composedApplicative = Applicative[IO] compose Applicative[ErrorOr]

  implicit def workflow: Case.Aux[Workflow, Map[String, Cwl] => Parse[Cwl]] =
    at[Workflow] {
      workflow =>
        initialMap: Map[String, Cwl] =>

        //Gather up all the filenames from the "run" section of workflow steps.
        val fileNames =
          workflow.
            steps.
            toList.
            flatMap(_.run.select[String].toList).
            filterNot(initialMap.keys.toList.contains)
          
        //read the files, parse them, and put them in a Map
        val cwlList: ParseValidated[List[(String, Cwl)]] =
          fileNames.traverse[ParseValidated, (String, Cwl)](CwlDecoder.decodeCwlAsValidated)


        //we have a list of tuples, what we actually want is a map
        val cwlMap: IO[ErrorOr[Map[String, Cwl]]] =
          //use the functor inherent in the Applicative to map over cwlList
          composedApplicative.map(cwlList)(_.toMap)

        EitherT { cwlMap.map(_.toEither) }.map{
          fileNameToCwl =>

            //Replace each step that links to another CWL file with the CWL instance
            val newWorkflow = lens[Workflow].steps.modify(workflow){
              _.map{ step =>
                lens[WorkflowStep].run.modify(step)(_.fold(RunToEmbeddedCwl).apply(fileNameToCwl ++ initialMap))
              }
            }

            newWorkflow.asCwl
        }
    }

  implicit def commandLineTool: Case.Aux[CommandLineTool, Map[String, Cwl] => Parse[Cwl]] =
    at[CommandLineTool] {
      clt => 
        Function.const(Monad[Parse].pure(clt.asCwl))
    }

  implicit def expressionTool: Case.Aux[ExpressionTool, Map[String, Cwl] => Parse[Cwl]] =
    at[ExpressionTool] {
      et =>
        Function.const(Monad[Parse].pure(et.asCwl))
    }
}
