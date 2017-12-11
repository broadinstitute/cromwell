package cwl

import cats.data.EitherT
import common.Checked
import cwl.CwlDecoder.Parse
import cwl.command.ParentName
import shapeless.Poly1
import common.validation.Checked._

private [cwl] object FlattenCwlFile extends Poly1 {
  private implicit val parentName = ParentName.empty
  
  private object CwlIdFromCwl extends Poly1 {
    implicit def fromWorkflow: Case.Aux[Workflow, String] = at[Workflow] { _.id }
    implicit def fromCommandLineTool: Case.Aux[CommandLineTool, String] = at[CommandLineTool] { _.id }
    implicit def fromExpressionLineTool: Case.Aux[ExpressionTool, String] = at[ExpressionTool] { _.id.get }
  }
  
  private def findRoot(cwls: Array[Cwl], root: String): Checked[Cwl] = {
    cwls
      .find(cwl => FileAndId(cwl.fold(CwlIdFromCwl)).id == root) match {
      case Some(rootCwl) => rootCwl.validNelCheck
      case None => s"Cannot find root object with id $root".invalidNelCheck
    }
  }
  
  private [cwl] object CwlFileRoot extends Poly1 {
    implicit def fromCwl: Case.Aux[Cwl, Option[String] => Checked[Cwl]] = at[Cwl] { cwl => _ => cwl.validNelCheck }
    implicit def fromCwlArray: Case.Aux[Array[Cwl], Option[String] => Checked[Cwl]] = at[Array[Cwl]] { cwls => 
      _ match {
        case Some(root) => findRoot(cwls, root)
        case None => "The supplied file contains an array of Cwl objects but no root has been provided. Please provided the cwl root".invalidNelCheck
      }
    }
  }

  type FileNameAndRootId = (String, Option[String])
  
  implicit def fromCwl: Case.Aux[Cwl, FileNameAndRootId => Parse[Cwl]] = at[Cwl] { cwl: Cwl => _: FileNameAndRootId => cwl.fold(AddEmbeddedCwl).apply(Map.empty[String, Cwl]) }
  
  implicit def fromCwlArray: Case.Aux[Array[Cwl], FileNameAndRootId => Parse[Cwl]] = at[Array[Cwl]] { cwls: Array[Cwl] => 
    fileAndRoot: FileNameAndRootId =>
      def processRoot(root: String): Parse[Cwl] = {
        // Regroup all parsed Cwls from this file in a map 
        val knownCwls: Map[WorkflowStepInputId, Cwl] = cwls.map(cwl => cwl.fold(CwlIdFromCwl) -> cwl).toMap
        // Find the root and use the map as a set of initially known cwls before embedding objects in the root
        findRoot(cwls, root) match {
          case Right(cwlRoot) => cwlRoot.fold(AddEmbeddedCwl).apply(knownCwls)
          case Left(error) => EitherT.leftT(error)
        }
      }
      
      fileAndRoot match {
        case (_, Some(root)) => processRoot(root)
        case _ => throw new Exception("The submission contained an array of CWL entities but no root was specified. Please specify the root object to be run when submitting multiple objects in the same file.")
      }
  }
}
