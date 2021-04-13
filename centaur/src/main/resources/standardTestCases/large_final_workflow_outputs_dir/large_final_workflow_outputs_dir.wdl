version 1.0

workflow large_final_workflow_outputs_dir {
    output {
        # Regression test CROM-6708 by batch copying a large file between two buckets.
        # In this case we're copying by using final_workflow_outputs_dir functionality.
        #
        # Because the file used in the test is large, via the workflow options we copy to
        # gs://cloud-cromwell-dev-self-cleaning-fast which is setup with a short lifecycle for deletion of objects.
        #
        # See also https://github.com/broadinstitute/rawls/blob/c39049945867d9d6d1bb5e1cbda30a09a19147f7/automation/src/test/scala/org/broadinstitute/dsde/test/api/RawlsApiSpec.scala#L768-L783
        #
        # HapMap is an old and very public project: https://en.wikipedia.org/wiki/International_HapMap_Project
        #
        # The link below are using the GCP Cloud Life Sciences public datasets:
        # https://cloud.google.com/life-sciences/docs/resources/public-datasets
        #
        # If someone at GCP pulls this copy of these reference files, there should be other public copies.
        # One example is the Broad's copy: gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz
        #
        # Also if/when we switch to multi-cloud, this file should be publicly available on other clouds, such as
        # https://registry.opendata.aws/broad-references/
        File out = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
        String patience = "It will take a few minutes while cromwell monitors the GCS copy of the workflow outputs"
    }

    # Final outputs are only copied if a call is executed creating a `backendAssignments`.
    # https://github.com/broadinstitute/cromwell/blob/58/engine/src/main/scala/cromwell/engine/workflow/lifecycle/finalization/CopyWorkflowOutputsActor.scala#L90
    call just_need_a_call
}

task just_need_a_call {
    command <<< true >>>
    # For now running on local to make the test run faster.
    # One day we may have docker image caching and push vs. poll status to make Papi jobs as fast as Local.
    runtime { backend: "Local" }
}
