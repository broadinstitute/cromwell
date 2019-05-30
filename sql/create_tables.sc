interp.configureCompiler(_.settings.YpartialUnification.value = true)
interp.repositories() ++= Seq(coursier.maven.MavenRepository(
  "https://oss.sonatype.org/content/repositories/releases"))


@

import $plugin.$ivy.`org.spire-math::kind-projector:0.9.7`
import scala.concurrent.duration._

/*

Migration of (Meta)Data from Cromwell's SQL datastore to Google Datastore


Algorithm:

every N seconds,
lastRowPulled <- query datastore by auto Id desc to get last row pushed to the Datastore
newData <- query Cromwell metadata from lastRow until end of all metadata
entities <- convert data into Datastore objects known as "Entities"
push the entites into Datastore

possible optimizations to be made:
Only pull data from datastore if we see there is no more data to be pulled from MySQL
pull data from mysql and load to DS in separate threads
*/

import $ivy.`org.tpolecat::doobie-core:0.7.0`
import $ivy.`mysql:mysql-connector-java:8.0.11`
import $ivy.`com.chuusai::shapeless:2.3.3`
import $ivy.`org.hsqldb:hsqldb:2.4.1`
/*
import $ivy.`org.scodec::scodec-core:1.11.3`
import $ivy.`org.scodec::scodec-bits:1.1.11`
import $ivy.`org.typelevel::cats-effect:1.3.1`
*/

import shapeless.syntax.std.tuple._

import doobie._
import doobie.implicits._
import cats.effect.IO
import cats.effect.IO._

import fs2.Stream
import java.util.concurrent.Executors
import java.time.Instant
import scala.{Stream => _}
import java.math.BigInteger

import cats.effect.Effect
import cats.instances.long._
import cats.instances.vector._
import cats.data.State
import cats.effect._
import cats.implicits._

import java.util.concurrent.Executors
import scala.concurrent.ExecutionContext

val executor = Executors.newCachedThreadPool()
val executionContext = ExecutionContext.fromExecutor(executor)
implicit val cs = IO.contextShift(executionContext)
implicit val timer = IO.timer(executionContext)

def blockingThreadPool[F[_]](implicit F: Sync[F]): Resource[F, ExecutionContext] =
  Resource(F.delay {
    val executor = Executors.newCachedThreadPool()
    val ec = ExecutionContext.fromExecutor(executor)
    (ec, F.delay(executor.shutdown()))
  })

val mysql = "jdbc:mysql://localhost/DatabaseName?rewriteBatchedStatements=true&useSSL=false"
val mem = "jdbc:hsqldb:mem:mymemdb"
val file = "jdbc:hsqldb:file:metadata;shutdown=false;hsqldb.tx=mvcc"

val db = mysql

val xa = Transactor.fromDriverManager[IO]("org.hsqldb.jdbcDriver", db, "ChooseAName","YourOtherPassword")


val create = sql"""
CREATE TABLE workflow (
    key BINARY(16),
    languageVersion_fk BINARY(16),
    backend_fk BINARY(16),
    name VARCHAR(65535),
    root VARCHAR(65535),
    submission DATETIME
);

CREATE TABLE run (
  workflow_fk BINARY(16),
  start DATETIME
);

CREATE TABLE workflow_input (
    key BINARY(16),
    workflow_fk BINARY(16),
    call_fk BINARY(16),
    name VARCHAR(65535)
)

CREATE TABLE task (
    key BINARY(16),
    call_caching_mode ENUM(''),
    docker_tag BINARY(16),
    fail_on_stderr BOOLEAN
)

CREATE TABLE backend {
  name VARCHAR(4096),
  key BINARY(16)
)

INSERT INTO task values(1,'',1,true);

""".update

val x = xa.yolo
import x._

create.run.void.transact(xa).attempt.unsafeRunSync

val putStrLn:String => ConnectionIO[Unit] = (s:String) => FC.delay(println(s))

val read = sql"SELECT call_caching_mode from task".query[String]
val y : Stream[ConnectionIO, Unit] = read.stream.evalMap(name => putStrLn(name))
y.compile.drain.transact(xa).attempt.unsafeRunSync
