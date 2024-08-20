package cromwell.database.slick.tables

import cromwell.database.sql.tables.GroupMetricsEntry
import slick.lifted.ProvenShape

import java.sql.Timestamp

trait GroupMetricsEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class GroupMetricsEntries(tag: Tag) extends Table[GroupMetricsEntry](tag, "GROUP_METRICS_ENTRY") {

    def groupMetricsEntryId = column[Long]("GROUP_METRICS_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def groupId = column[String]("GROUP_ID", O.Length(255))

    def quotaExhaustionDetected = column[Timestamp]("QUOTA_EXHAUSTION_DETECTED")

//    def baseProjection = (groupId, quotaExhaustionDetected)

    /*
       In Slick, adding '?' at end of column lifts it into an Option column. The reason the auto-incremented primary key
       in Slick is represented as Option[T], for example 'groupMetricsEntryId' in this case, is because when a row is
       inserted into the table the value for the primary key is unknown (since it's auto-incremented by DB) and hence
       in this case the value for that column is represented as 'None'. After insertion, database assigns it an
       auto-incremented value and at this point the value for that column is represented as 'Some(value)'.
     */
    override def * : ProvenShape[GroupMetricsEntry] = (groupId,
                                                       quotaExhaustionDetected,
                                                       groupMetricsEntryId.?
    ) <> ((GroupMetricsEntry.apply _).tupled, GroupMetricsEntry.unapply)

//    def forUpdate = baseProjection.shaped <> (
//      tuple => GroupMetricsEntry.tupled(tuple :+ None),
//      GroupMetricsEntry.unapply(_: GroupMetricsEntry).map(_.reverse.tail.reverse)
//    )

    def ixGroupMetricsEntryGi = index("IX_GROUP_METRICS_ENTRY_GI", groupId, unique = false)
  }

  protected val groupMetricsEntries = TableQuery[GroupMetricsEntries]

  val groupMetricsEntryIdsAutoInc = groupMetricsEntries returning groupMetricsEntries.map(_.groupMetricsEntryId)

  val quotaExhaustionForGroupId = Compiled((groupId: Rep[String]) =>
    for {
      groupMetricsEntry <- groupMetricsEntries
      if groupMetricsEntry.groupId === groupId
    } yield groupMetricsEntry.quotaExhaustionDetected
  )
}

// CompiledFunction[driver.api.Rep[String] => Query[Rep[String], String, Seq], driver.api.Rep[String], String, Query[Rep[String], String, Seq], Seq[String]]
// CompiledFunction[driver.api.Rep[String] => Query[Rep[String], String, Seq], driver.api.Rep[String], String, Query[Rep[String], String, Seq], Seq[String]]