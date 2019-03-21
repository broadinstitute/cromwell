package example

import java.io.{InputStream, OutputStream}

import cats.syntax.either._
import com.google.gson.Gson
import com.github.gekomad.ittocsv.core.ToCsv._
import _root_.io.circe.generic.auto._
import _root_.io.circe._
import _root_.io.circe.parser._
import _root_.io.circe.syntax._
import _root_.io.circe.parser.{decode => cDecode}
import java.math.BigInteger
import java.nio.charset.StandardCharsets
import java.sql.{PreparedStatement, ResultSet, Timestamp}
import java.time.Instant
import java.{lang, util}

import org.apache.beam.sdk.io.jdbc.JdbcIO
import cats.kernel.Monoid
import cats.data.State
import com.github.gekomad.ittocsv.core.{CsvStringEncoder, ParseFailure}
import com.google.gson.internal.LinkedTreeMap
//import com.spotify.scio._
import org.apache.beam.sdk
import org.apache.beam.sdk.io.{FileIO, FileSystems, TextIO}
import org.apache.beam.sdk.options.{PipelineOptions, PipelineOptionsFactory}
import org.apache.beam.sdk.Pipeline
import org.apache.beam.sdk.coders.{StringUtf8Coder, StructuredCoder}
import org.apache.beam.sdk.io.FileIO.{MatchAll, MatchConfiguration, ReadMatches}
import org.apache.beam.sdk.io.fs.MatchResult.Metadata
import org.apache.beam.sdk.io.fs.{EmptyMatchTreatment, MatchResult}
import org.apache.beam.sdk.io.jdbc.JdbcIO.PreparedStatementSetter
import org.apache.beam.sdk.transforms.Combine.CombineFn
import org.apache.beam.sdk.transforms.DoFn.ProcessElement
import org.apache.beam.sdk.transforms.windowing.BoundedWindow
import org.apache.beam.sdk.transforms._
import org.apache.beam.sdk.values.ValueInSingleWindow.Coder
import org.apache.beam.sdk.values._
import org.apache.beam.vendor.guava.v20_0.com.google.common.io.ByteStreams

import scala.collection.JavaConverters._

/*
sbt "runMain [PACKAGE].WordCount
  --project=[PROJECT] --runner=DataflowRunner --zone=[ZONE]
  --input=gs://dataflow-samples/shakespeare/kinglear.txt
  --output=gs://[BUCKET]/[PATH]/wordcount"
*/

object WordCount {
  import Csv._
  implicit val csvFormat = com.github.gekomad.ittocsv.parser.IttoCSVFormat.default

  def func[A,B](f: A => B): ProcessFunction[A,B] = new ProcessFunction[A,B] {
    override def apply(input: A): B = f(input)
  }

  def main(cmdlineArgs: Array[String]): Unit = {
    //val host = "localhost"
    val host = "google"

    //TODO Pas
    case class GoogleDB(user: String, pass: String, name: String, zone: String, project: String)
    def dbString(db: GoogleDB): String = s"jdbc:mysql://$host/cromwell?cloudSqlInstance=${db.project}:${db.zone}:${db.name}&socketFactory=com.google.cloud.sql.mysql.SocketFactory&user=${db.user}&password=${db.pass}&useUnicode=true&characterEncoding=UTF-8"
    val anicholsClone = GoogleDB("","","","", "")
    val mcovarrClone = GoogleDB("","","","", "")
    val metadata3 = GoogleDB("","","","", "")
    val papiV2 = GoogleDB("","","","", "")

    val readDb = papiV2
    val writeDb = metadata3

    val jdbcP = Pipeline.create(PipelineOptionsFactory.fromArgs(cmdlineArgs:_*).create())

    //val end = 250000L
    val start = 4L
    val end = 150L

    val  inputList: util.List[BigInteger] = List.range[Long](start, end).map(java.math.BigInteger.valueOf).map((java.math.BigInteger.valueOf(10000)).multiply).asJava

    val mysqlDriverString = "com.mysql.jdbc.Driver"

    jdbcP.apply(Create.of(inputList))
    .apply(JdbcIO.readAll().
        withDataSourceConfiguration(JdbcIO.DataSourceConfiguration.create(mysqlDriverString,dbString(readDb))).
        withQuery("SELECT * FROM METADATA_ENTRY WHERE METADATA_JOURNAL_ID > ? LIMIT 10000").
        withParameterSetter(new PreparedStatementSetter[BigInteger] {
          override def setParameters(element: BigInteger, preparedStatement: PreparedStatement): Unit = preparedStatement.setObject(1, element)
        }).
        withCoder(StringUtf8Coder.of()).
        withRowMapper(new JdbcIO.RowMapper[String] {
          override def mapRow(resultSet: ResultSet): String = {
            import resultSet._
            val tuple = (
              getObject(1, classOf[java.math.BigInteger]),
              getString(2),
              getString(3),
              getString(4),
              getInt(5),
              getInt(6),
              getString(7),
              getTimestamp(8).toString,
              getString(9)
            )
            //toCsv(tuple)
            val gson = new Gson()
            gson.toJson(tuple)
          }})).
      apply(JdbcIO.write().
        withDataSourceConfiguration(JdbcIO.DataSourceConfiguration.create(mysqlDriverString,dbString(writeDb))).
        withStatement("insert into METADATA_ENTRY values(?,?,?,?,?,?,?,?,?)").
        withBatchSize(10000).
        withPreparedStatementSetter(new PreparedStatementSetter[String] {
          override def setParameters(element: String, preparedStatement: PreparedStatement): Unit =
            element match {
              case string =>
                import com.github.gekomad.ittocsv.core.FromCsv._

                //This is a sad necessity of using Gson, the only usable serialization library I can get to work
                //Maybe should try Protobuf next
                val gson = new Gson()
                val result = gson.fromJson(string, classOf[MetadataColumns])
                import preparedStatement._
                val v5Double:Double  =  result._5.asInstanceOf[java.lang.Double]
                val v5 = v5Double.toInt
                val v6Double: Double = result._6.asInstanceOf[java.lang.Double]
                val v6 = v6Double.toInt
                val v8String = result._8
                val v8 = Timestamp.valueOf(v8String)

                setObject(1,result._1)
                    setString(2,result._2)
                    setString(3,result._3)
                    setString(4,result._4)
                    setInt(5,v5)
                    setInt(6,v6)
                    setString(7,result._7)
                    setTimestamp(8,v8)
                    setString(9,result._9)
            }}))
    jdbcP.run().waitUntilFinish()
  }
}

case class MetaToString() extends DoFn[Metadata, String] {
  @ProcessElement
  def processElement(c: ProcessContext) {
  val write: Metadata = c.element();
  val size = write.toString
  c.output(size)
}
}

abstract class StateCombiner[I, S, O](implicit monoid: Monoid[S], state: Option[I] => State[S, O])extends CombineFn[I, S, O] {

  def createAccumulator(): S = monoid.empty

  def addInput(accumulator: S, input: I): S = state(Option(input)).run(accumulator).value._1

  def mergeAccumulators(accumulators: lang.Iterable[S]): S = monoid.combineAll(accumulators.asScala)

  def extractOutput(accumulator: S): O = state(None).run(accumulator).value._2
}

case class CirceCoder[A]()(implicit encoder: Encoder[A], decoder: Decoder[A]) extends StructuredCoder[A] {

  import org.apache.beam.sdk.util.VarInt;
  override def encode(value: A, outStream: OutputStream): Unit = {
    Option(value).foreach {value =>
      val e: Json = value.asJson
      val s = Printer.noSpaces.pretty(e)
      val bytes = s.getBytes
      VarInt.encode(bytes.size, outStream)
      outStream.write(bytes)
    }
  }

  override def decode(inStream: InputStream): A = {
    val numBytes = VarInt.decodeInt(inStream)
    val bytes = new scala.Array[Byte](numBytes)
    ByteStreams.readFully(inStream,bytes)
    val json = new String(bytes, StandardCharsets.UTF_8);
    cDecode[A](json).getOrElse(throw new RuntimeException("unable to parse"))
  }

  override def getCoderArguments: util.List[_ <: sdk.coders.Coder[_]] = new java.util.ArrayList()

  override def verifyDeterministic(): Unit = ()
}

object Csv {
  implicit def bigIntegerEncoder: CsvStringEncoder[java.math.BigInteger] = createEncoder(_.toString)
  implicit def bigInteger: String => Either[ParseFailure, java.math.BigInteger] = s => Either.catchNonFatal(new BigInteger(s)).leftMap(t => ParseFailure(t.getMessage))

  implicit def instantEncoder: CsvStringEncoder[Instant] = createEncoder(_.toString)
  implicit def instant: String => Either[ParseFailure, Instant] = s => Either.catchNonFatal(Instant.parse(s)).leftMap(t => ParseFailure(t.getMessage))
}


