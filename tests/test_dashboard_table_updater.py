import pytest

from genie import dashboard_table_updater as dash_update


@pytest.mark.parametrize(
    "input_string_time, expected_output_time",
    [
        ("2018-10-25T20:16:07.959Z", 1540498567000),
        ("2018-10-25T20:16:07.959", 1540498567000),
        ("2018-04-06T18:30:00", 1523039400000),
    ],
    ids=["utc_time", "local_time_zone", "no_milliseconds"],
)
def test_that_string_to_unix_epoch_time_milliseconds_gives_expected_time(
    input_string_time, expected_output_time
):
    output = dash_update.string_to_unix_epoch_time_milliseconds(input_string_time)
    assert output == expected_output_time
