package cromwell.util

object TerminalUtil {
  def highlight(colorCode:Int, string:String) = s"\033[38;5;${colorCode}m$string\033[0m"
  def mdTable(rows: Seq[Seq[String]], header: Seq[String]): String = {
    def maxWidth(lengths: Seq[Seq[Int]], column: Int) = lengths.map { length => length(column) }.max
    val widths = (rows :+ header).map { row => row.map { s => s.length } }
    val maxWidths = widths.head.indices.map { column => maxWidth(widths, column) }
    val tableHeader = header.indices.map { i => header(i).padTo(maxWidths(i), " ").mkString("") }.mkString("|")
    val tableDivider = header.indices.map { i => "-" * maxWidths(i) }.mkString("|")
    val tableRows = rows.map { row =>
      val mdRow = row.indices.map { i => row(i).padTo(maxWidths(i), " ").mkString("") }.mkString("|")
      s"|$mdRow|"
    }
    s"|$tableHeader|\n|$tableDivider|\n${tableRows.mkString("\n")}\n"
  }
}
