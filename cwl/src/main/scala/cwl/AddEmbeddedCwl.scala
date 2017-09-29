package cwl

import shapeless.{Path => _, _}
import cats.syntax.traverse._
import cats.instances.list._
import cats.effect.IO
import CwlDecoder.Parse
import cats.Monad
import cats.data.EitherT._
import cats.data.EitherT
import lenthall.validation.ErrorOr._

object AddEmbeddedCwl extends Poly1 {

  import CwlDecoder._

  implicit def workflow =
    at[Workflow] {
      workflow =>

        //Gather up all the filenames from the "run" section of workflow steps.
        val fileNames =
          workflow.
            steps.
            toList.
            flatMap(_.run.select[String].toList)

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
                lens[WorkflowStep].run.modify(step)(_.fold(RunToEmbeddedCwl).apply(fileNameToCwl))
              }
            }

            newWorkflow.asCwl
        }
    }

  implicit def commandLineTool =
    at[CommandLineTool] {
      clt =>
        Monad[Parse].pure(clt.asCwl)
    }
}
