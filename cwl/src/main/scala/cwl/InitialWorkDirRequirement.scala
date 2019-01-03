package cwl

import cwl.InitialWorkDirRequirement._
import eu.timepit.refined.W
import shapeless.{:+:, CNil, _}

final case class InitialWorkDirRequirement(
                                            `class`: W.`"InitialWorkDirRequirement"`.T,
                                            listing: IwdrListing
                                          ) {
  val listings: Array[IwdrListingArrayEntry] = listing.fold(IwdrListingArrayPoly)

  override def toString: WorkflowStepInputId =
    s"""InitialWorkDirRequirement(
      |  ${listings.mkString(System.lineSeparator + "  ")}
      |)""".stripMargin
}

/**
  *  Short for "Directory Entry"
  *  @see <a href="http://www.commonwl.org/v1.0/CommandLineTool.html#Dirent">Dirent Specification</a>
  *
  *  Split into two cases because entryName is only optional if entry is an Expression
  */
trait Dirent {
  def writable: Option[Boolean]
  def writableWithDefault = writable.getOrElse(false)
}

final case class ExpressionDirent(
                                   entry: Expression,
                                   entryname: Option[StringOrExpression],
                                   writable: Option[Boolean]
                                 ) extends Dirent

final case class StringDirent(
                               entry: String,
                               entryname: StringOrExpression,
                               writable: Option[Boolean]
                             ) extends Dirent

object InitialWorkDirRequirement {

  // "ExpressionDirent" has to come before StringDirent because expressions are matched by "String" first if we do it the
  // other way round.
  final type IwdrListingArrayEntry = File :+: Directory :+: ExpressionDirent :+: StringDirent :+: StringOrExpression :+: CNil
  final type IwdrListing = Array[IwdrListingArrayEntry] :+: StringOrExpression :+: CNil

  object IwdrListingArrayPoly extends Poly1 {
    implicit val caseArrayIwdrListingArrayEntry: Case.Aux[Array[IwdrListingArrayEntry], Array[IwdrListingArrayEntry]] = {
      at {
        identity
      }
    }

    implicit val caseStringOrExpression: Case.Aux[StringOrExpression, Array[IwdrListingArrayEntry]] = {
      at {
        stringOrExpression =>
          Array(Coproduct[IwdrListingArrayEntry](stringOrExpression))
      }
    }
  }
}
