from datetime import timedelta, datetime, timezone
import dateutil.parser
from numpy import exp, log, polyfit
import os
import requests


def main(call_endpoint = True):
    time_series_response = call_time_series_endpoint() if call_endpoint else read_response_from_file()
    points = time_series_response["timeSeries"][0]["points"]
    data = [datum_from_raw_time_point(p) for p in points]

    # Sorting is not strictly necessary but convenient when sanity checking, and I am pro-sanity.
    data.sort(key = lambda d: d.timestamp)

    # Fitting
    #
    # y = a * e ^ (b * x)
    #
    # can be like fitting
    #
    # log y = log a + b * x
    #
    # if logarithms are taken of the y values

    # Months and TiB seem the most natural units to describe the growth of this db.
    # 30 days in a month as far as this script is concerned.
    seconds_per_month = (60 * 60 * 24 * 30)
    bytes_per_tib = (1 << 40)

    now_utc = datetime.now(timezone.utc)
    x = [(d.timestamp - now_utc).total_seconds() / seconds_per_month for d in data]
    y = [log(d.bytes / bytes_per_tib) for d in data]

    # numpy's polyfit is perverse, returns coefficients in *descending* order of degree.
    [b, log_a] = polyfit(x, y, 1)
    a = exp(log_a)
    print("For <db size in TiB> = a * exp(b * <months from now>), a and b are estimated to be {a} and {b} respectively.".format(**locals()))
    size_from_now = a * exp(b * 12)
    print("The Cromwell database size is predicted to be {size_from_now:.2f} TiB 12 months from now.".format(**locals()))

    maximum_database_size_tib = 30

    # y = <db size in TiB>, x = <months from now>
    # To find out when the database will explode, solve for x and set y to the maximum size.
    # y = a * e^(b*x)          # original
    # ln(y) = ln(a) + b*x      # take logarithms of both sides
    # ln(y) - ln(a) = b*x      # subtract ln(b) from both sides
    # (ln(y) - ln(a))/b = x    # divide both sides by c
    # x = (ln(y) - ln(a))/b    # swap sides
    when_explode_months = (log(maximum_database_size_tib) - log(a)) / b
    print("Given a maximum size of {maximum_database_size_tib:.2f} TiB, it is estimated that the Cromwell database will explode in {when_explode_months:.2f} months. Please plan accordingly.".format(**locals()))


def datum_from_raw_time_point(time_point):
    string_timestamp = time_point["interval"]["startTime"]
    bytes = time_point["value"]["doubleValue"]
    return Datum(dateutil.parser.parse(string_timestamp), bytes)


def read_response_from_file():
    import json
    line = open("response.json", "r").readlines()[0].strip()
    time_series_response = json.loads(line)
    return time_series_response


def call_time_series_endpoint():
    google_project = os.getenv("CROMWELL_GOOGLE_PROJECT")
    cloudsql_instance = os.getenv("CROMWELL_CLOUDSQL_INSTANCE")
    access_token = os.getenv("CROMWELL_METRICS_ACCESS_TOKEN")

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
    return response.json()


def predict_fullness():
    maximum_database_size_tib = os.getenv("CROMWELL_MAXIMUM_DATABASE_SIZE_TIB")
    database_forecast_limit_months = os.getenv("CROMWELL_DATABASE_FORECAST_LIMIT_MONTHS")


class Datum:
    def __init__(self, timestamp, bytes):
        self.timestamp = timestamp
        self.bytes = bytes


if __name__ == '__main__':
    # print("Hello from Python")
    # print(os.environ)
    main(call_endpoint = False)
