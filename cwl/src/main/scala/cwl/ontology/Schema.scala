package cwl.ontology

import java.util.concurrent.Executors

import cats.effect.IO
import cats.syntax.traverse._
import cats.instances.list._
import com.google.common.cache.{Cache, CacheBuilder}
import com.typesafe.config.{Config, ConfigFactory}
import com.typesafe.scalalogging.Logger
import common.util.IORetry
import common.util.IORetry.StatefulIoError
import common.validation.ErrorOr._
import common.validation.Validation._
import cwl.ontology.Schema._
import mouse.all._
import net.ceedubs.ficus.Ficus._
import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.model._
import org.semanticweb.owlapi.model.parameters.OntologyCopy
import org.semanticweb.owlapi.reasoner.structural.StructuralReasonerFactory
import org.semanticweb.owlapi.reasoner.{OWLReasoner, OWLReasonerFactory}
import org.semanticweb.owlapi.util.OWLAPIStreamUtils
import org.slf4j.LoggerFactory

import scala.collection.JavaConverters._
import scala.concurrent.ExecutionContext
import scala.util.Try

/**
  * OWL/RDF Schema lookup.
  *
  * @param schemaIris IRI paths to OWL/RDF schemas.
  * @param namespaces Additional/Override namespace prefixes.
  */
case class Schema(schemaIris: Seq[String],
                  namespaces: Map[String, String],
                  ontologyManager: OWLOntologyManager = OWLManager.createOWLOntologyManager,
                  reasonerFactory: OWLReasonerFactory = new StructuralReasonerFactory) {

  /**
    * Returns a full IRI based on a full or abbreviated IRI.
    *
    * A full IRI wrapped in < and > will not be checked for abbreviations but will be returned without the wrapper.
    *
    * @see https://www.w3.org/TR/owl2-syntax/#IRIs
    */
  def fullIri(name: String): String = getIri(name).getIRIString

  /**
    * Returns true if child is an equivalent of ancestor, and if not then recursively checks if any of the child's
    * super classes are an equivalent of ancestor.
    */
  def isSubClass(child: String, ancestor: String): Boolean = {
    val childIri: IRI = getIri(child)
    val ancestorIri: IRI = getIri(ancestor)
    val childClass: OWLClass = dataFactory.getOWLClass(childIri)
    val ancestorClass: OWLClass = dataFactory.getOWLClass(ancestorIri)
    val schemaReasoner: OWLReasoner = reasonerFactory.createReasoner(schemaOntology)
    try {
      Schema.isSubClass(schemaReasoner, childClass, ancestorClass)
    } finally {
      schemaReasoner.dispose()
    }
  }

  private val dataFactory: OWLDataFactory = ontologyManager.getOWLDataFactory
  private val schemaOntology: OWLOntology = ontologyManager.createOntology()
  private val schemaPrefixManager: PrefixManager =
    ontologyManager.getOntologyFormat(schemaOntology).asPrefixOWLDocumentFormat

  {
    def addToSchema(originalOntology: OWLOntology): Unit = {
      schemaOntology.addAxioms(originalOntology.axioms())
      val originalOntologyFormat: OWLDocumentFormat = ontologyManager.getOntologyFormat(originalOntology)
      if (originalOntologyFormat.isPrefixOWLDocumentFormat)
        schemaPrefixManager.copyPrefixesFrom(originalOntologyFormat.asPrefixOWLDocumentFormat)
    }

    val errorOr: ErrorOr[Unit] = for {
      ontologies <- schemaIris.toList.traverse(loadOntologyFromIri(ontologyManager))
      _ = ontologies.foreach(addToSchema)
    } yield ()
    errorOr.toTry("Error loading schemas").get

    // Add any namespace overrides
    namespaces foreach {
      case (prefixName, prefix) => schemaPrefixManager.setPrefix(prefixName, prefix)
    }
  }

  /**
    * Returns the full IRI for the name using the prefixManager.
    *
    * TODO: Not 100% sure why you can't ask owl-api for an iri-with-prefix and have it looked up automatically.
    *
    * There does seem to be a difference between abbreviated and full IRIs in the official spec, where full IRIs are
    * wrapped in < >. But this doesn't seem to be the format used by CWL nor owl-api.
    *
    * Please update this comment if/when one knows the correct behavior.
    *
    * @see https://www.w3.org/TR/owl2-syntax/#IRIs
    */
  private def getIri(name: String): IRI = {
    Try(schemaPrefixManager.getIRI(name)).getOrElse(IRI.create(name))
  }
}

object Schema {
  // Extending StrictLogging creates a circular dependency here for some reason, so making the logger ourselves
  private val logger: Logger = Logger(LoggerFactory.getLogger(getClass.getName))
  private [ontology] val ontologyConfig = ConfigFactory.load.as[Config]("ontology")
  private val ontologyConfiguration = OntologyConfiguration(ontologyConfig)
  private [ontology] val cacheConfig = ontologyConfig.getAs[Config]("cache")

  // Simple cache to avoid reloading the same ontologies too often
  private val ontologyCache = cacheConfig.map(makeOntologyCache)

  private implicit val statefulIoError = StatefulIoError.noop[Unit]
  private implicit val timer = cats.effect.IO.timer(ExecutionContext.fromExecutor(Executors.newFixedThreadPool(ontologyConfiguration.poolSize)))

  private [ontology] def makeOntologyCache(config: Config): Cache[IRI, OWLOntology] = {
    val cacheConfig = CacheConfiguration(config)
    logger.info(s"Ontology cache size: ${cacheConfig.maxSize}")
    CacheBuilder.newBuilder()
      .maximumSize(cacheConfig.maxSize)
      .build[IRI, OWLOntology]()
  }

  /**
    * Returns the absolute path for a file, possibly relative to parent.
    */
  def getIriPath(parent: String, path: String): String = IRI.create(parent).resolve(path).getIRIString

  /**
    * Load an ontology either from an IRI.
    */
  private [ontology] def loadOntologyFromIri(ontologyManager: OWLOntologyManager, cache: Option[Cache[IRI, OWLOntology]] = ontologyCache)(schemaIri: String): ErrorOr[OWLOntology] = {
    validate {
      val iri = IRI.create(schemaIri)
      cache.flatMap(_.getIfPresent(iri) |> Option.apply) match {
        case Some(ontology) =>
          ontologyManager.copyOntology(ontology, OntologyCopy.DEEP)
        case _ =>
          logger.info(s"Loading ${iri.toURI.toString}")
          val ontology = loadOntologyFromIri(ontologyManager, iri)
          cache.foreach(_.put(iri, ontology))
          ontology
      }
    }
  }

  // Loading the ontology can fail transiently, so put retires around it. See https://github.com/protegeproject/webprotege/issues/298
  private [ontology] def loadOntologyFromIri(ontologyManager: OWLOntologyManager, iri: IRI): OWLOntology = {
    val load = IO { ontologyManager.loadOntologyFromOntologyDocument(iri) }
    IORetry.withRetry[OWLOntology, Unit](load, (), ontologyConfiguration.retries, ontologyConfiguration.backoff).unsafeRunSync()
  }

  /**
    * Returns true if child is an equivalent of ancestor, and if not then recursively checks if any of the child's
    * super classes are an equivalent of ancestor.
    */
  private def isSubClass(reasoner: OWLReasoner, childClass: OWLClass, ancestorClass: OWLClass): Boolean = {
    val equivalent: Set[OWLClass] = reasoner.getEquivalentClasses(childClass).asScala.toSet + childClass
    if (equivalent.contains(ancestorClass)) {
      true
    } else {
      val parentClasses: Set[OWLClass] = for {
        equivalentClass <- equivalent
        parentClass <- OWLAPIStreamUtils
          .asSet(reasoner.getSuperClasses(equivalentClass).entities)
          .asScala
          .toSet[OWLClass]
      } yield parentClass
      parentClasses.collect({
        case superClass: OWLClass if isSubClass(reasoner, superClass, ancestorClass) => superClass
      }).nonEmpty
    }
  }
}
