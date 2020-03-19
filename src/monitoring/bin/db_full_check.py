from datetime import timedelta, datetime, timezone
import dateutil.parser
import numpy
import os
import requests


def main():
    time_series_response = call_time_series_endpoint()
    points = time_series_response.json()["timeSeries"][0]["points"]
    time_series_data = [Datum(dateutil.parser.parse(p["interval"]["startTime"]), p["value"]["doubleValue"]) for p in points]

    # Sorting is not strictly necessary but convenient when sanity checking.
    time_series_data.sort(key = lambda d: d.timestamp)

    for datum in time_series_data:
        print(datum.timestamp, datum.bytes)


def call_time_series_endpoint():
    google_project = os.getenv("CROMWELL_GOOGLE_PROJECT")
    cloudsql_instance = os.getenv("CROMWELL_CLOUDSQL_INSTANCE")
    access_token = os.getenv("METRICS_ACCESS_TOKEN")

    now_utc = datetime.now(timezone.utc)
    now = datetime.isoformat(now_utc)
    then = datetime.isoformat(now_utc - timedelta(days = 60))

    # print("now is '%s' and then is '%s'"%(now, then))

    url='https://monitoring.clients6.google.com/v3/projects/{google_project}/timeSeries'.format(**locals())
    filter='metric.type="cloudsql.googleapis.com/database/disk/bytes_used" AND (resource.label.database_id="{google_project}:{cloudsql_instance}")'.format(**locals())

    # print("url is", url, "and filter is", filter)
    response = requests.get(
        url,
        params = {
            'filter': filter,
            'interval.startTime': then,
            'interval.endTime': now,
            'aggregation.alignmentPeriod': '+10800s',
            'aggregation.perSeriesAligner': 'ALIGN_MEAN'
        },
        headers = {
            'Authorization': 'Bearer ' + access_token
        }
    )
    return response


def predict_fullness():
    maximum_database_size_tib = os.getenv("CROMWELL_MAXIMUM_DATABASE_SIZE_TIB")
    database_forecast_limit_months = os.getenv("CROMWELL_DATABASE_FORECAST_LIMIT_MONTHS")


class Datum:
    def __init__(self, timestamp, bytes):
        self.timestamp = timestamp
        self.bytes = bytes

if __name__ == '__main__':
    print("Hello from Python")
    print(os.environ)
    main()
