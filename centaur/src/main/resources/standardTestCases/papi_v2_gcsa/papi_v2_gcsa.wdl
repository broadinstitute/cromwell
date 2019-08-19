version 1.0

workflow papi_v2_gcsa {
    call get_token_info
    output {
        String email = get_token_info.json.email
        String scope = get_token_info.json.scope
    }
}

# Confirm that even though the service account (SA) specified by the backend configuration creates the pipeline job,
# instead the google compute service account (GCSA) workflow option is the actual account used to run the container.
# https://cloud.google.com/genomics/reference/rest/Shared.Types/Metadata#VirtualMachine.FIELDS.service_account
# https://cromwell.readthedocs.io/en/stable/wf_options/Google/#google-pipelines-api-workflow-options
task get_token_info {
    command <<<
        apt-get install --assume-yes jq > /dev/null
        curl \
            --fail \
            --silent \
            --request POST \
            "https://www.googleapis.com/oauth2/v2/tokeninfo?access_token=$(gcloud auth application-default print-access-token)" \
        | jq \
            --monochrome-output \
            '{email, scope: .scope | split(" ") | sort | join(" ")}'
    >>>

    runtime {
        docker: "google/cloud-sdk:slim"
        backend: "papi-v2-gcsa"
    }

    output {
        Object json = read_json(stdout())
    }
}
