# Based on https://stackoverflow.com/a/39895650
from __future__ import print_function
from docker import Client
import os

cli = Client(base_url='unix://var/run/docker.sock')

our_hostname = os.environ.get("HOSTNAME")

all_containers = cli.containers()

# filter out ourself by hostname
our_container = \
  [c for c in all_containers if c['Id'][:12] == our_hostname[:12]][0]

project_name = our_container['Labels']['com.docker.compose.project']
service_name = 'cromwell'
filters = [
  'com.docker.compose.project={}'.format(project_name)
]

containers = cli.containers(filters={'label': filters})

hostname_list = []
for container in [x for x in containers if x['Labels']['com.docker.compose.service'].startswith(service_name)]:
  project = container['Labels']['com.docker.compose.project']
  service = container['Labels']['com.docker.compose.service']
  number = container['Labels']['com.docker.compose.container-number']
  hostname = '_'.join([project, service, number])
  hostname_list.append(hostname)
  print(hostname)
# end for
