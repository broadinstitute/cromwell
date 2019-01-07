#!/usr/bin/env python

from functools import reduce
from google.cloud.monitoring_v3 import MetricServiceClient
from google.cloud.monitoring_v3.types import LabelDescriptor, MetricDescriptor, TimeSeries
from os import environ
import psutil as ps
import requests
from signal import signal, SIGTERM
from sys import stderr
from time import sleep, time

def get_metadata(key):
  return requests.get(
    'http://metadata.google.internal/computeMetadata/v1/instance/' + key,
    headers={'Metadata-Flavor': 'Google'}
  ).text

def reset():
  global memory_used, disk_used, disk_reads, disk_writes, report_time

  # Explicitly reset the CPU counter, because the first call of this method always reports 0
  ps.cpu_percent()

  memory_used = 0

  disk_used = 0
  disk_reads = disk_io('read_count')
  disk_writes = disk_io('write_count')

  report_time = 0

def measure():
  global memory_used, disk_used, report_time

  memory_used = max(memory_used, MEMORY_SIZE - mem_usage('available'))
  disk_used = max(disk_used, disk_usage('used'))

  report_time += MEASUREMENT_TIME_SEC
  sleep(MEASUREMENT_TIME_SEC)

def mem_usage(param):
  return getattr(ps.virtual_memory(), param)

def disk_usage(param):
  return reduce(
    lambda usage, mount: usage + getattr(ps.disk_usage(mount), param),
    DISK_MOUNTS, 0,
  )

def disk_io(param):
  return getattr(ps.disk_io_counters(), param)

def format_gb(value_bytes):
  return '%.1f' % round(value_bytes / 2**30, 1)

def get_metric(key, value_type, unit, description):
  return client.create_metric_descriptor(PROJECT_NAME, MetricDescriptor(
    type='/'.join(['custom.googleapis.com', METRIC_ROOT, key]),
    description=description,
    metric_kind='GAUGE',
    value_type=value_type,
    unit=unit,
    labels=LABEL_DESCRIPTORS,
  ))

def create_time_series(series):
  client.create_time_series(PROJECT_NAME, series)

def get_time_series(metric_descriptor, value):
  series = TimeSeries()

  series.metric.type = metric_descriptor.type
  labels = series.metric.labels
  labels['workflow_id'] = WORKFLOW_ID
  labels['task_call_name'] = TASK_CALL_NAME
  labels['task_call_index'] = TASK_CALL_INDEX
  labels['task_call_attempt'] = TASK_CALL_ATTEMPT
  labels['cpu_count'] = CPU_COUNT_LABEL
  labels['mem_size'] = MEMORY_SIZE_LABEL
  labels['disk_size'] = DISK_SIZE_LABEL

  series.resource.type = 'gce_instance'
  series.resource.labels['zone'] = ZONE
  series.resource.labels['instance_id'] = INSTANCE

  point = series.points.add(value=value)
  point.interval.end_time.seconds = int(time())

  return series

def report():
  create_time_series([
    get_time_series(CPU_UTILIZATION_METRIC, { 'double_value': ps.cpu_percent() }),
    get_time_series(MEMORY_UTILIZATION_METRIC, { 'double_value': memory_used / MEMORY_SIZE * 100 }),
    get_time_series(DISK_UTILIZATION_METRIC, { 'double_value': disk_used / DISK_SIZE * 100 }),
    get_time_series(DISK_READS_METRIC, { 'double_value': (disk_io('read_count') - disk_reads) / report_time }),
    get_time_series(DISK_WRITES_METRIC, { 'double_value': (disk_io('write_count') - disk_writes) / report_time }),
  ])

### Define constants

# Cromwell variables passed to the container
# through environmental variables
WORKFLOW_ID = environ['WORKFLOW_ID']
TASK_CALL_NAME = environ['TASK_CALL_NAME']
TASK_CALL_INDEX = environ['TASK_CALL_INDEX']
TASK_CALL_ATTEMPT = environ['TASK_CALL_ATTEMPT']
DISK_MOUNTS = environ['DISK_MOUNTS'].split()

# GCP instance name, zone and project
# from instance introspection API
INSTANCE = get_metadata('name')
_, PROJECT, _, ZONE = get_metadata('zone').split('/')

client = MetricServiceClient()
PROJECT_NAME = client.project_path(PROJECT)

METRIC_ROOT = 'wdl_task'

MEASUREMENT_TIME_SEC = 1
REPORT_TIME_SEC = 60

LABEL_DESCRIPTORS = [
  LabelDescriptor(
    key='workflow_id',
    description='Cromwell workflow ID',
  ),
  LabelDescriptor(
    key='task_call_name',
    description='Cromwell task call name',
  ),
  LabelDescriptor(
    key='task_call_index',
    description='Cromwell task call index',
  ),
  LabelDescriptor(
    key='task_call_attempt',
    description='Cromwell task call attempt',
  ),
  LabelDescriptor(
    key='cpu_count',
    description='Number of virtual cores',
  ),
  LabelDescriptor(
    key='mem_size',
    description='Total memory size, GB',
  ),
  LabelDescriptor(
    key='disk_size',
    description='Total disk size, GB',
  ),
]

CPU_COUNT = ps.cpu_count()
CPU_COUNT_LABEL = str(CPU_COUNT)

MEMORY_SIZE = mem_usage('total')
MEMORY_SIZE_LABEL = format_gb(MEMORY_SIZE)

DISK_SIZE = disk_usage('total')
DISK_SIZE_LABEL = format_gb(DISK_SIZE)

CPU_UTILIZATION_METRIC = get_metric(
  'cpu_utilization', 'DOUBLE', '%',
  '% of CPU utilized in a Cromwell task call',
)

MEMORY_UTILIZATION_METRIC = get_metric(
  'mem_utilization', 'DOUBLE', '%',
  '% of memory utilized in a Cromwell task call',
)

DISK_UTILIZATION_METRIC = get_metric(
  'disk_utilization', 'DOUBLE', '%',
  '% of disk utilized in a Cromwell task call',
)

DISK_READS_METRIC = get_metric(
  'disk_reads', 'DOUBLE', '{reads}/s',
  'Disk read IOPS in a Cromwell task call',
)

DISK_WRITES_METRIC = get_metric(
  'disk_writes', 'DOUBLE', '{writes}/s',
  'Disk write IOPS in a Cromwell task call',
)

### Detect container termination

def signal_handler(signum, frame):
  global running
  running = False

running = True
signal(SIGTERM, signal_handler)

### Main loop
#
# It continuously measures runtime metrics every MEASUREMENT_TIME_SEC,
# and reports them to Stackdriver Monitoring API every REPORT_TIME_SEC.
#
# However, if it detects a container termination signal,
# it *should* report the final metric
# right after the current measurement, and then exit normally.

reset()
while running:
  measure()
  if not running or report_time >= REPORT_TIME_SEC:
    report()
    reset()
exit(0)
