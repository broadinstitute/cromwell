package cromwell.database.sql.joins

import cromwell.database.sql.tables._

case class CallCachingJoin
(
  callCachingEntry: CallCachingEntry,
  callCachingHashEntries: Seq[CallCachingHashEntry],
  callCachingAggregationEntry: Option[CallCachingAggregationEntry],
  callCachingSimpletonEntries: Seq[CallCachingSimpletonEntry],
  callCachingDetritusEntries: Seq[CallCachingDetritusEntry]
)
