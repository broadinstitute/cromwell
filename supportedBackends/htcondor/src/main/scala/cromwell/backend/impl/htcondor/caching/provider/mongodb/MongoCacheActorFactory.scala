package cromwell.backend.impl.htcondor.caching.provider.mongodb

import akka.actor.Props
import com.mongodb.casbah.{MongoClient, MongoCollection}
import com.typesafe.config.Config
import cromwell.backend.impl.htcondor.caching.CacheActorFactory

class MongoCacheActorFactory(config: Config) extends CacheActorFactory {
  val dbHost = config.getString("cache.db.host")
  val dbPort = config.getInt("cache.db.port")
  val dbName = config.getString("cache.db.name")
  val dbCollectionName = config.getString("cache.db.collection")
  val dbInstance: MongoClient = MongoClient(dbHost, dbPort)
  val db = dbInstance(dbName)
  val collection: MongoCollection = db(dbCollectionName)

  override def getCacheActorProps(forceRewrite: Boolean): Props = {
    Props(new MongoCacheActor(collection, forceRewrite))
  }

}
