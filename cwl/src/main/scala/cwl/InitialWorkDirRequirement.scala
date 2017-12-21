package cwl

import cwl.InitialWorkDirRequirement._
import eu.timepit.refined.W
import shapeless.{:+:, CNil, _}

final case class InitialWorkDirRequirement(
                                            `class`: W.`"InitialWorkDirRequirement"`.T,
                                            listing: IwdrListing
                                          ) {
  val listings: Array[IwdrListingArrayEntry] = listing match {
    case IwdrListingArray(array) => array
    case StringOrExpression(soe) => Array(Coproduct[IwdrListingArrayEntry](soe))
  }

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

  final type IwdrListingArrayEntry = File :+: Directory :+: ExpressionDirent :+: StringDirent :+: StringOrExpression :+: CNil
  final type IwdrListing = Array[IwdrListingArrayEntry] :+: StringOrExpression :+: CNil

  object IwdrListingArrayEntry {
    object StringDirent {
      def unapply(e: IwdrListingArrayEntry): Option[(String, StringOrExpression, Boolean)] =
        e.select[StringDirent].map(sd => (sd.entry, sd.entryname, sd.writableWithDefault))
    }
    object ExpressionDirent {
      def unapply(e: IwdrListingArrayEntry): Option[(Expression, Option[StringOrExpression], Boolean)] =
        e.select[ExpressionDirent].map(sd => (sd.entry, sd.entryname, sd.writableWithDefault))
    }
  }

  object IwdrListingArray {
    def unapply(listing: IwdrListing): Option[Array[IwdrListingArrayEntry]] =
      listing.select[Array[IwdrListingArrayEntry]]
  }
}
