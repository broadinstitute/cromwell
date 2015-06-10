package cromwell.util

object TerminalUtil {
  def highlight(colorCode:Int, string:String) = s"\033[38;5;${colorCode}m${string}\033[0m"
}
