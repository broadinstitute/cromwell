package cromwell.engine

import cromwell.binding._
import cromwell.binding.values._
import org.apache.commons.codec.digest.DigestUtils

import scala.collection.immutable.TreeMap


trait Hashable extends Any {
  def md5Sum: String
}

object Hashing {

  implicit class HashableString(val value: String) extends AnyVal with Hashable {
    def md5Sum: String = DigestUtils.md5Hex(value)
  }

  implicit class HashableWdlValue(val wdlValue: WdlValue) extends AnyVal {

    private def symbolHash(hash: String) = SymbolHash((wdlValue.getClass.getCanonicalName + hash).md5Sum)

    private def symbolHash[K](hashedMap: Map[K, SymbolHash])(implicit ord: Ordering[K]): SymbolHash = {
      // productIterator returns an Iterator over the elements of a Tuple2 Map entry.
      val concatenatedMap = TreeMap(hashedMap.toArray: _*) flatMap { _.productIterator } mkString ""
      symbolHash(concatenatedMap)
    }

    def getHash(descriptor: WorkflowDescriptor): Option[SymbolHash] =
      if (descriptor.configCallCaching) Option(computeHash(descriptor.fileHasher)) else None

    def computeHash(implicit hasher: FileHasher): SymbolHash = {
      wdlValue match {
        case w: WdlObject => symbolHash(w.value mapValues { _.computeHash(hasher) })
        case w: WdlMap => symbolHash(w.value map { case (k, v) => k.computeHash(hasher) -> v.computeHash(hasher) })
        case w: WdlArray => symbolHash(w.value map { _.computeHash(hasher) } mkString "")
        case w: WdlFile => hasher(w)
        case w => symbolHash(w.valueString)
      }
    }
  }

}
