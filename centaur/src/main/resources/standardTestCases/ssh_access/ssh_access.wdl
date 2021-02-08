version 1.0

workflow ssh_access {
    call ssh

    output {
        String out = ssh.out
    }
}

task ssh {
    meta {
        volatile: true
    }
    command <<<
        set -euo pipefail

        CURL_ARGS=(--silent --show-error --fail --header "Metadata-Flavor: Google")

        # Retrieve the host, zone, and project from the metadata
        GCE_HOST="$(curl "${CURL_ARGS[@]}" http://metadata.google.internal/computeMetadata/v1/instance/name)"
        GCE_ZONE="$(curl "${CURL_ARGS[@]}" http://metadata.google.internal/computeMetadata/v1/instance/zone | cut -d/ -f4)"
        GCE_PROJECT="$(curl "${CURL_ARGS[@]}" http://metadata.google.internal/computeMetadata/v1/project/project-id)"

        # Loopback doesn't seem to work unless we allow access from our own IP address.
        # Tried:
        #   - gcloud compute ssh --internal-ip
        #     - Ironically, ssh using internal ips *does* work between two different instances, but not on loopback.
        #   - gcloud compute firewall-rules create --source-service-accounts=[our SA]
        #     - Doesn't work targeting all IP addresses or targeting only [our SA]
        #   - gcloud compute ssh --dry-run & replacing the ip with a local iP
        #     - 127.0.0.1 is the IP of the docker container
        #     - but this container needs to connect over to the other ssh-server container
        #     - trying to grab the `ip addr eth0` also seemed to be different across two different docker containers

        # Instead let's temporarily modify the firewall allowing access:
        #   - from the IP address we're holding
        #   - to any machines started by our service account
        #   - that are also listening on port 22

        # Retrieve our external IP address
        GCE_IP="$(
            gcloud \
                --quiet \
                compute instances describe \
                --zone "$GCE_ZONE" \
                --project "$GCE_PROJECT" \
                "$GCE_HOST" \
                --format='get(networkInterfaces[0].accessConfigs[0].natIP)'
        )"

        # Retrieve the account we're logged in as assuming it was the SA used to create this machine
        GCE_ACCOUNT="$(gcloud config get-value account)"

        # Tag the name of the firewall rule with our ip address. If we're getting collisions then our cleanup is broken!
        GCE_FIREWALL_RULE="centaur-temp-${GCE_IP//./-}"

        # Create the temporary firewall rule
        gcloud \
            --quiet \
            compute firewall-rules create \
            "${GCE_FIREWALL_RULE}" \
            --description="created $(date +%s) for centaur test ssh_access.test" \
            --allow=tcp:22 \
            --priority=1022 \
            --source-ranges="${GCE_IP}"/32 \
            --target-service-accounts="${GCE_ACCOUNT}"

        # Create a trap to delete the temporary firewall rule
        delete_firewall() {
            gcloud \
                --quiet \
                compute firewall-rules delete \
                "${GCE_FIREWALL_RULE}"
        }
        trap delete_firewall EXIT

        # Run the actual test.
        # If someone understands why the ssh host keys are sometimes other-than-expected please explain in a comment and
        # possibly remove the `--strict-host-key=no`.
        gcloud \
            --quiet \
            compute ssh \
            --zone "$GCE_ZONE" \
            --project "$GCE_PROJECT" \
            "$GCE_HOST" \
            --strict-host-key-checking=no \
            --command="echo tt1276104 > $PWD/over_ssh.txt"
    >>>

    output {
        String out = read_string("over_ssh.txt")
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk"
    }
}
