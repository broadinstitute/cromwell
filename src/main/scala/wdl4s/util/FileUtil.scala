package wdl4s.util

import scala.util.{Failure, Success, Try}

object FileUtil {
  // FIXME: Keep?
  def parseTsv(tsv: String): Try[Array[Array[String]]] = {
    val table = tsv.split("\n").map(_.split("\t"))
    table.map(_.length).toSet match {
      case s if s.size > 1 => Failure(new UnsupportedOperationException("TSV is not uniform"))
      case _ => Success(table)
    }
  }
}
