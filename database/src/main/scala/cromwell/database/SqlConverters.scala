package cromwell.database

import java.sql.{Clob, Timestamp}
import java.util.Date
import javax.sql.rowset.serial.SerialClob

object SqlConverters {

  implicit class DateToTimestamp(val date: Date) extends AnyVal {
    def toTimestamp: Timestamp = date match {
      case timestamp: Timestamp => timestamp
      case _ => new Timestamp(date.getTime)
    }
  }

  implicit class ClobToRawString(val clob: Clob) extends AnyVal {
    def toRawString: String = clob.getSubString(1, clob.length.toInt) // yes, it starts at 1
  }

  implicit class StringToClob(val str: String) extends AnyVal {
    def toClob: Clob = new SerialClob(str.toCharArray)

    def toNonEmptyClob: Option[Clob] = if (str.isEmpty) None else Option(new SerialClob(str.toCharArray))
  }

}
