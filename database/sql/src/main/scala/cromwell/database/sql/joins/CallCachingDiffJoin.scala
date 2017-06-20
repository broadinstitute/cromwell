package cromwell.database.sql.joins

import cromwell.database.sql.tables.CallCachingEntry

case class CallCachingDiffJoin(cacheEntryA: CallCachingEntry, cacheEntryB: CallCachingEntry, diff: Seq[(Option[(String, String)], Option[(String, String)])])
