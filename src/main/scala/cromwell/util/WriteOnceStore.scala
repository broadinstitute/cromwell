package cromwell.util

import scala.collection.concurrent.TrieMap
import scala.util.{Success, Failure, Try}

class WriteOnceStore[K, V] {
  private val writeOnceStore = new TrieMap[K, V]()

  /**
   * Insert a key value pair into the store. If the key already exists, returns a Failure, otherwise it will return
   * a Success containing an immutable copy of the current store w/ the key inserted
   */
  def insert(key: K, value: V): Try[Map[K, V]] = {
    println(s"$key => $value")
    writeOnceStore.putIfAbsent(key, value) match {
      case Some(x) => Failure(new IllegalStateException(s"$key was already present in the store"))
      case None => Success(toMap)
    }
  }

  /**
   * @return An immutable copy of the current store map
   */
  def toMap: Map[K, V] = writeOnceStore.readOnlySnapshot().toMap
}