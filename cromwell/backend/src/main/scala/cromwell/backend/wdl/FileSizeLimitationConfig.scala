package cromwell.backend.wdl

import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.readers.ValueReader
import net.ceedubs.ficus.Ficus._

trait FileSizeLimitationConfig {

  def readLinesLimit: Int

  def readBoolLimit: Int

  def readIntLimit: Int

  def readFloatLimit: Int

  def readStringLimit: Int

  def readJsonLimit: Int

  def readTsvLimit: Int

  def readMapLimit: Int

  def readObjectLimit: Int
}

object FileSizeLimitationConfig {
  private val config = ConfigFactory.load.getConfig("system")

  def fileSizeLimitationConfig: FileSizeLimitationConfig = config.as[FileSizeLimitationConfig]("input-read-limits")

  implicit val configReader : ValueReader[FileSizeLimitationConfig] = ValueReader.relative{c =>
    def f(s: String) = c.as[Int](s)
    new FileSizeLimitationConfig {
      val readLinesLimit =  f("lines")
      val readBoolLimit =   f("bool")
      val readIntLimit =    f("int")
      val readFloatLimit =  f("float")
      val readStringLimit = f("string")
      val readJsonLimit =   f("json")
      val readTsvLimit =    f("tsv")
      val readMapLimit =    f("map")
      val readObjectLimit = f("object")
    }
  }
}

