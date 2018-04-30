package cwl.ontology

import cats.instances.list._
import cats.syntax.traverse._
import common.validation.ErrorOr._
import common.validation.Validation._
import cwl.ontology.Schema._
import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.model._
import org.semanticweb.owlapi.reasoner.structural.StructuralReasonerFactory
import org.semanticweb.owlapi.reasoner.{OWLReasoner, OWLReasonerFactory}
import org.semanticweb.owlapi.util.OWLAPIStreamUtils

import scala.collection.JavaConverters._
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
      ontologies <- schemaIris.toList.traverse[ErrorOr, OWLOntology](loadOntologyFromIri(ontologyManager))
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

  /**
    * Returns the absolute path for a file, possibly relative to parent.
    */
  def getIriPath(parent: String, path: String): String = IRI.create(parent).resolve(path).getIRIString

  /**
    * Load an ontology either from an IRI.
    */
  private def loadOntologyFromIri(ontologyManager: OWLOntologyManager)(schemaIri: String): ErrorOr[OWLOntology] = {
    validate {
      val iri = IRI.create(schemaIri)
      ontologyManager.loadOntologyFromOntologyDocument(iri)
    }
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
