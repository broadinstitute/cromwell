package wom.callable

import wom.types.{WomMapType, WomStringType}
import wom.values.{WomMap, WomString, WomValue}

/**
  * Parameter documentation quoted from CWL Spec.
  *
  * @param outputPath runtime.outdir: an absolute path to the designated output directory
  * @param tempPath runtime.tmpdir: an absolute path to the designated temporary directory
  * @param cores runtime.cores: number of CPU cores reserved for the tool process
  * @param ram runtime.ram: amount of RAM in mebibytes (2**20) reserved for the tool process
  * @param outputPathSize runtime.outdirSize: reserved storage space available in the designated output directory
  * @param tempPathSize runtime.tmpdirSize: reserved storage space available in the designated temporary directory
  */
case class RuntimeEnvironment(outputPath: String,
                              tempPath: String,
                              cores: Int,
                              ram: Double,
                              outputPathSize: Long,
                              tempPathSize: Long) {

  def cwlMap: WomValue = {

    val womMap: Map[WomValue, WomValue] = Map(
      "outdir" -> outputPath,
      "tmpdir" -> tempPath,
      "cores" -> cores.toString,
      "ram" -> ram.toString,
      "outdirSize" -> outputPathSize.toString,
      "tmpdirSize" -> tempPathSize.toString
    ).map{
      case (key, value) => WomString(key) -> WomString(value)
    }

    WomMap(
      WomMapType(WomStringType, WomStringType),
      womMap
    )
  }

}

