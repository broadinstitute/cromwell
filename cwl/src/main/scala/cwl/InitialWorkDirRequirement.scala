package cwl

import cwl.InitialWorkDirRequirement.{ExpressionOrString, IwdrListingArrayEntry, _}
import eu.timepit.refined.W
import shapeless.{:+:, CNil, _}

final case class InitialWorkDirRequirement(
                                            `class`: W.`"InitialWorkDirRequirement"`.T,
                                            listing: Array[IwdrListingArrayEntry] :+: ExpressionOrString :+: CNil
                                          ) {
  val listings: Array[IwdrListingArrayEntry] = listing match {
    case IwdrListingArray(array) => array
    case ExpressionOrString(eos) => Array(Coproduct[IwdrListingArrayEntry](eos))
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
                                   entryname: Option[ExpressionOrString],
                                   writable: Option[Boolean]
                                 ) extends Dirent

final case class StringDirent(
                               entry: String,
                               entryname: ExpressionOrString,
                               writable: Option[Boolean]
                             ) extends Dirent

object InitialWorkDirRequirement {

  final type IwdrListingArrayEntry = File :+: Directory :+: StringDirent :+: ExpressionDirent :+: ExpressionOrString :+: CNil

  object IwdrListingArrayEntry {
    object StringDirent {
      def unapply(e: IwdrListingArrayEntry): Option[(String, ExpressionOrString, Boolean)] =
        e.select[StringDirent].map(sd => (sd.entry, sd.entryname, sd.writableWithDefault))
    }
  }

  object IwdrListingArray {
    def unapply(listing: Array[IwdrListingArrayEntry] :+: ExpressionOrString :+: CNil): Option[Array[IwdrListingArrayEntry]] =
      listing.select[Array[IwdrListingArrayEntry]]
  }
  object Expression { def unapply(listing: Array[IwdrListingArrayEntry] :+: ExpressionOrString :+: CNil): Option[Array[IwdrListingArrayEntry]] = listing.select[Array[IwdrListingArrayEntry]] }



  final type ExpressionOrString = Expression :+: String :+: CNil
  object ExpressionOrString {
    def unapply(listing: Array[IwdrListingArrayEntry] :+: ExpressionOrString :+: CNil): Option[ExpressionOrString] = listing.select[ExpressionOrString]
    object Expression { def unapply(eos: ExpressionOrString): Option[Expression] = eos.select[Expression] }
    object String { def unapply(eos: ExpressionOrString): Option[String] = eos.select[String] }
  }
}
