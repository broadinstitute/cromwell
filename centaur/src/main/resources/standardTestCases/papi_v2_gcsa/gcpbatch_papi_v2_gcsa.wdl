version 1.0

workflow papi_v2_gcsa {
    call get_token_info
    output {
        String email = get_token_info.json.email
        String scopes = get_token_info.json.scopes
        File service_account = get_token_info.service_account
    }
}

# Confirm that even though the service account (SA) specified by the backend configuration creates the pipeline job,
# instead the google compute service account (GCSA) workflow option is the actual account used to run the container.
# https://cloud.google.com/genomics/reference/rest/Shared.Types/Metadata#VirtualMachine.FIELDS.service_account
# https://cromwell.readthedocs.io/en/stable/wf_options/Google/#google-pipelines-api-workflow-options
task get_token_info {
    command <<<
        apt-get install --assume-yes jq > /dev/null

        curl --fail --silent \
            'http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/?recursive=true' \
            -H 'Metadata-Flavor: Google' > service_accounts.json

        cat service_accounts.json | jq --monochrome-output \
            '.default | {email, scopes: .scopes | sort | join(" ")}'
    >>>

    runtime {
        docker: "mirror.gcr.io/google/cloud-sdk:slim"
        backend: "GCPBATCH-gcsa"
    }

    output {
        Object json = read_json(stdout())
        File service_account = "service_accounts.json"
    }
}
