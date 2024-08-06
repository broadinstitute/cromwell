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

    override def * : ProvenShape[GroupMetricsEntry] = (groupId,
                                                       quotaExhaustionDetected,
                                                       groupMetricsEntryId.?
    ) <> ((GroupMetricsEntry.apply _).tupled, GroupMetricsEntry.unapply)

    def ixGroupMetricsEntryGi = index("IX_GROUP_METRICS_ENTRY_GI", groupId, unique = false)
  }

  protected val groupMetricsEntries = TableQuery[GroupMetricsEntries]

  val groupMetricsEntryIdsAutoInc = groupMetricsEntries returning groupMetricsEntries.map(_.groupMetricsEntryId)

}
