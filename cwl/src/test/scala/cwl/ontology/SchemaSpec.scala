package cwl.ontology

import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import org.semanticweb.owlapi.apibinding.OWLManager

class SchemaSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "Schema"

  val resolveTests = Table(
    ("base", "iri", "expected"),
    ("file:///hello/world", "custom.owl", "file:/hello/custom.owl"),
    ("file:///hello/world", "file:./custom.owl", "file:./custom.owl"),
    ("file:///hello/world", "file:///custom.owl", "file:///custom.owl"),
    ("http://hello/world", "custom.owl", "http://hello/custom.owl"),
    ("http://hello/world", "file:./custom.owl", "file:./custom.owl"),
    ("http://hello/world", "file:///custom.owl", "file:///custom.owl"),
  )

  forAll(resolveTests) { (base, iri, expected) =>
    it should s"resolve $iri relative to $base to create $expected" in {
      Schema.getIriPath(base, iri) should be(expected)
    }
  }

  val fullIriTests = Table(
    ("schemaDescription", "schema", "iri", "fullIri"),
    ("galaxy", galaxySchema, "edam:format_1929", "http://edamontology.org/format_1929"),
    ("galaxy", galaxySchema, "<edam:format_1929>", "edam:format_1929"),
    ("galaxy", galaxySchema, "edam:format_2200", "http://edamontology.org/format_2200"),
    ("galaxy", galaxySchema, "http://galaxyproject.org/formats/fasta", "http://galaxyproject.org/formats/fasta"),
    ("galaxy", galaxySchema, "gx:fasta", "http://galaxyproject.org/formats/fasta"),
    ("galaxy", galaxySchema, "not:real", "not:real"),
    ("galaxy", galaxySchema, "owl:Thing", "http://www.w3.org/2002/07/owl#Thing"),
    ("galaxy", galaxySchema, "owl:Nothing", "http://www.w3.org/2002/07/owl#Nothing"),
    ("edam", edamSchema, "edam:format_1929", "http://purl.obolibrary.org/obo/edam#format_1929"),
    ("edam w/ namespace", edamNamespaceSchema, "edam:format_1929", "http://edamontology.org/format_1929"),
  )

  forAll(fullIriTests) { (schemaDescription, schema, iri, fullIri) =>
    it should s"expand $iri to $fullIri using schema $schemaDescription" in {
      schema.fullIri(iri) should be(fullIri)
    }
  }

  val subClassTests = Table(
    ("isSubClass", "schemaDescription", "schema", "child", "ancestor"),
    (true, "galaxy", galaxySchema, "edam:format_1929", "edam:format_1929"),
    (true, "galaxy", galaxySchema, "edam:format_1929", "http://galaxyproject.org/formats/fasta"),
    (true, "galaxy", galaxySchema, "http://galaxyproject.org/formats/fasta", "edam:format_1929"),
    (true, "galaxy", galaxySchema, "gx:fasta", "edam:format_1929"),
    (true, "galaxy", galaxySchema, "gx:fasta", "edam:format_2330"),
    (false, "galaxy", galaxySchema, "not:real", "edam:format_1929"),
    (false, "galaxy", galaxySchema, "edam:format_1929", "not:real"),
    (false, "galaxy", galaxySchema, "owl:Thing", "owl:Nothing"),
    (true, "galaxy", galaxySchema, "owl:Nothing", "owl:Thing"),
    (false, "galaxy", galaxySchema, "owl:Thing", "edam:format_1929"),
    (true, "galaxy", galaxySchema, "edam:format_1929", "owl:Thing"),
    (false, "galaxy", galaxySchema, "edam:format_1929", "owl:Nothing"),
    (true, "galaxy", galaxySchema, "owl:Nothing", "edam:format_1929"),
    (true, "galaxy", galaxySchema, "edam:format_1929", "http://edamontology.org/format_2330"),
    (false, "edam", edamSchema, "edam:format_1929", "http://edamontology.org/format_2330"),
    (true, "edam w/ namespace", edamNamespaceSchema, "edam:format_1929", "http://edamontology.org/format_2330"),
  )

  forAll(subClassTests) { (isSubClass, schemaDescription, schema, child, ancestor) =>
    it should s"return $isSubClass that $child is a sub class of $ancestor when checking schema $schemaDescription" in {
      schema.isSubClass(child, ancestor) should be(isSubClass)
    }
  }
  
  it should "cache ontologies" in {
    val cache = Schema.makeOntologyCache(ConfigFactory.parseString("max-size=5"))
    val ontologyManager1 = OWLManager.createOWLOntologyManager
    val ontologyManager2 = OWLManager.createOWLOntologyManager
    val ontology1 = Schema.loadOntologyFromIri(ontologyManager1, Option(cache))(getClass.getResource("EDAM.owl").toExternalForm)
    ontology1.isValid shouldBe true
    val ontology2 = Schema.loadOntologyFromIri(ontologyManager2, Option(cache))(getClass.getResource("EDAM.owl").toExternalForm)
    ontology2.isValid shouldBe true
    ontology2.getOrElse(fail()).getOWLOntologyManager shouldBe ontologyManager2
  }

  private lazy val galaxySchema = {
    Schema(
      List(
        getClass.getResource("EDAM.owl").toExternalForm,
        getClass.getResource("gx_edam.ttl").toExternalForm
      ),
      Map.empty
    )
  }

  private lazy val edamSchema = {
    Schema(
      List(
        getClass.getResource("EDAM.owl").toExternalForm
      ),
      Map.empty
    )
  }

  private lazy val edamNamespaceSchema = {
    Schema(
      List(
        getClass.getResource("EDAM.owl").toExternalForm
      ),
      Map("edam" -> "http://edamontology.org/")
    )
  }
}
