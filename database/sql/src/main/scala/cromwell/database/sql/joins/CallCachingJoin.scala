package cromwell.database.sql.joins

import cromwell.database.sql.tables.{CallCachingDetritusEntry, CallCachingEntry, CallCachingHashEntry, CallCachingSimpletonEntry}

case class CallCachingJoin
(
  callCachingEntry: CallCachingEntry,
  callCachingHashEntries: Seq[CallCachingHashEntry],
  callCachingSimpletonEntries: Seq[CallCachingSimpletonEntry],
  callCachingDetritusEntries: Seq[CallCachingDetritusEntry]
)
