#!/usr/bin/env bash

cat << EOF >> /var/lib/postgresql/data/postgresql.conf
max_connections = 300
EOF
