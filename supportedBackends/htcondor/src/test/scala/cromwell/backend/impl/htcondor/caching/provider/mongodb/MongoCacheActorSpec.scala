package cromwell.backend.impl.htcondor.caching.provider.mongodb

import akka.actor.ActorSystem
import akka.testkit.{ImplicitSender, TestActorRef, TestKit}
import com.mongodb.casbah.MongoCollection
import com.mongodb.casbah.commons.MongoDBObject
import com.mongodb.util.JSON
import com.mongodb.{DBObject, WriteResult}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.{MemorySize, BackendJobDescriptorKey}
import cromwell.backend.BackendJobExecutionActor.JobSucceededResponse
import cromwell.backend.impl.htcondor.HtCondorRuntimeAttributes
import cromwell.backend.impl.htcondor.caching.CacheActor._
import cromwell.backend.impl.htcondor.caching.exception.CachedResultNotFoundException
import cromwell.backend.impl.htcondor.caching.provider.mongodb.model.{KryoSerializedObject, MongoCachedExecutionResult}
import cromwell.backend.impl.htcondor.caching.provider.mongodb.serialization.KryoSerDe
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.core.JobOutput
import org.mockito.Mockito
import org.mockito.Mockito._
import org.scalatest.mockito.MockitoSugar
import org.scalatest.{BeforeAndAfter, BeforeAndAfterAll, MustMatchers, WordSpecLike}
import wdl4s.TaskCall
import wdl4s.values.WdlString

class MongoCacheActorSpec extends TestKit(ActorSystem("MongoCacheProviderActorSpecSystem")) with WordSpecLike with MustMatchers
  with BeforeAndAfter with BeforeAndAfterAll with ImplicitSender with MockitoSugar with KryoSerDe {

  import spray.json._
  import cromwell.backend.impl.htcondor.caching.provider.mongodb.model.MongoCachedExecutionResultProtocol._

  val config: Config = ConfigFactory.load()
  val mongoDbCollectionMock = mock[MongoCollection]
  val memorySize = MemorySize.parse("0.512 GB").get
  val diskSize = MemorySize.parse("1.024 GB").get
  val runtimeConfig = HtCondorRuntimeAttributes(ContinueOnReturnCodeSet(Set(0)), Some("tool-name"), Some("/workingDir"), Some("/outputDir"), true, 1, memorySize, diskSize, None)
  val jobHash = "88dde49db10f1551299fb9937f313c10"
  val taskStatus = "done"
  val succeededResponseMock = JobSucceededResponse(BackendJobDescriptorKey(TaskCall(Option("taskName"), null, null, null), None, 0), None, Map("test" -> JobOutput(WdlString("Test"))), None, Seq.empty)
  val serSucceededRespMock = KryoSerializedObject(serialize(succeededResponseMock))
  val cachedExecutionResult = MongoCachedExecutionResult(jobHash, serSucceededRespMock)
  val cachedExecutionResultDbObject = JSON.parse(cachedExecutionResult.toJson.toString).asInstanceOf[DBObject]
  val query = MongoDBObject("hash" -> jobHash)

  after {
    Mockito.reset(mongoDbCollectionMock)
  }

  override def afterAll = shutdown()

  "A CacheActor" should {
    "return an ExecutionResultFound when read an execution result from cache" in {
      when(mongoDbCollectionMock.findOne(query)) thenReturn Some(cachedExecutionResultDbObject)
      val cacheActor = TestActorRef(new MongoCacheActor(mongoDbCollectionMock))
      cacheActor ! ReadExecutionResult(jobHash)
      expectMsg(ExecutionResultFound(succeededResponseMock))
      verify(mongoDbCollectionMock, atLeastOnce).findOne(query)
    }

    "return an ExecutionResultNotFound when it can't find an execution result in cache" in {
      when(mongoDbCollectionMock.findOne(query)) thenReturn None
      val cacheActor = TestActorRef(new MongoCacheActor(mongoDbCollectionMock))
      cacheActor ! ReadExecutionResult(jobHash)
      expectMsg(ExecutionResultNotFound)
      verify(mongoDbCollectionMock, atLeastOnce).findOne(query)
    }

    "return ExecutionResultStored when it stores an execution result" in {
      when(mongoDbCollectionMock.findOne(query)) thenThrow new CachedResultNotFoundException("")
      when(mongoDbCollectionMock.insert(cachedExecutionResultDbObject)) thenReturn new WriteResult(0, true, "")
      val cacheActor = TestActorRef(new MongoCacheActor(mongoDbCollectionMock))
      cacheActor ! StoreExecutionResult(jobHash, succeededResponseMock)
      expectMsg(ExecutionResultStored("88dde49db10f1551299fb9937f313c10"))
      verify(mongoDbCollectionMock, atLeastOnce).insert(cachedExecutionResultDbObject)
    }

    "return ExecutionResultAlreadyExist when it tries to store an existing execution result" in {
      when(mongoDbCollectionMock.findOne(query)) thenReturn Some(cachedExecutionResultDbObject)
      val cacheActor = TestActorRef(new MongoCacheActor(mongoDbCollectionMock))
      cacheActor ! StoreExecutionResult(jobHash, succeededResponseMock)
      expectMsg(ExecutionResultAlreadyExist)
      verify(mongoDbCollectionMock, atLeastOnce).findOne(query)
    }

    "return ExecutionResultNotFound when try to read and force re-write flag is enabled" in {
      when(mongoDbCollectionMock.findOne(query)) thenReturn Some(cachedExecutionResultDbObject)
      val cacheActor = TestActorRef(new MongoCacheActor(mongoDbCollectionMock, true))
      cacheActor ! ReadExecutionResult(jobHash)
      expectMsg(ExecutionResultNotFound)
      verify(mongoDbCollectionMock, atLeastOnce).findOne(query)
    }

    "return ExecutionResultStored when try to store and force re-write flag is enabled" in {
      when(mongoDbCollectionMock.findOne(query)) thenThrow new CachedResultNotFoundException("")
      when(mongoDbCollectionMock.insert(cachedExecutionResultDbObject)) thenReturn new WriteResult(0, true, "")
      when(mongoDbCollectionMock.remove(query)) thenReturn new WriteResult(0, true, "")
      val cacheActor = TestActorRef(new MongoCacheActor(mongoDbCollectionMock))
      cacheActor ! StoreExecutionResult(jobHash, succeededResponseMock)
      expectMsg(ExecutionResultStored("88dde49db10f1551299fb9937f313c10"))
      verify(mongoDbCollectionMock, atLeastOnce).insert(cachedExecutionResultDbObject)
      verify(mongoDbCollectionMock, atLeastOnce).remove(query)
    }
  }
}
