def calculate_overlap(i1_begin, i1_end, i2_begin, i2_end):
    """
    Calculates percentage of overlapping positions between two internvals, [i1_begin, i1_end[ and [i2_begin, i2_end[.

    @parameter i1_begin Begin position of interval 1.
    @parameter i1_end End position of interval 1. Should be greater or equal to i1_begin.
    @parameter i2_begin Begin position of interval 2.
    @parameter i2_end End position of interval 2. Should be greater or equal to i2_begin.

    @returns Percentage in a float number in the range [0.0, 100.0]
    """

    # Make sure the ranges are valid
    assert i1_begin <= i1_end
    assert i2_begin <= i2_end

    # Calculate the number of overlapping bases
    overlapping_bases = max(0, min(i1_end, i2_end) - max(i1_begin, i2_begin))

    # Calculate the size of the larger interval
    larger_interval_size = max(1, max(i1_end - i1_begin, i2_end - i2_begin))

    # Divide the two numbers number of overlapping bases with the size of larger interval
    return float(100.0 * overlapping_bases) / float(larger_interval_size)


# Unit tests for the function 'calculate_overlap()'
assert calculate_overlap(0, 1, 2, 3) == 0.0
assert calculate_overlap(2, 3, 0, 1) == 0.0
assert calculate_overlap(0, 1, 0, 1) == 100.0
assert calculate_overlap(0, 1, 1, 1) == 0.0
assert calculate_overlap(0, 2, 1, 1) == 0.0
assert calculate_overlap(0, 2, 1, 2) == 50.0
assert calculate_overlap(0, 3, 1, 2) == float(100)/float(3)
assert calculate_overlap(0, 100, 32, 62) == float(3000)/float(100)
assert calculate_overlap(0, 1000, 32, 62) == float(3000)/float(1000)
assert calculate_overlap(0, 1000, 32, 1042) == float(100.0 * (1000-32))/float(1010)
assert calculate_overlap(32, 1042, 0, 1000) == float(100.0 * (1000-32))/float(1010)


def calculate_mean(data):
    """
    Calculates the sample mean of data. Returns a float.
    """
    n = len(data)
    assert n > 0
    return float(sum(data)/float(len(data)))


# Unit tests for the function 'calculate_mean()'
assert calculate_mean([0.0]) == 0.0
assert calculate_mean([2]) == 2.0
assert calculate_mean([1, 2]) == 1.5
assert calculate_mean([100, -100]) == 0.0
assert calculate_mean(range(10)) == 4.5
assert calculate_mean(range(11)) == 5.0


def calculate_ssd(data):
    """
    Calculates the sum of square deviation of data.
    """
    c = calculate_mean(data)
    ss = sum((float(x) - c)**2 for x in data)
    return float(ss)


def calculate_stddev(data, degrees_of_freedom=1.0):
    """
    Calculates the population standard deviation of data.
    """
    n = len(data)

    if n < 2:
        return 0.0  # variance requires at least two data points

    ssd = calculate_ssd(data)
    pvar = ssd / float(n - degrees_of_freedom)
    return pvar ** 0.5


def calculate_stddev_no_outliers(data, degrees_of_freedom=1.0):
    """
    Calculates the population standard deviation of data but ignores outliers
    """
    n = len(data)

    if n < 2:
        return 0.0  # variance requires at least two data points

    c = calculate_mean(data)
    ssd_values = [(float(x) - c) ** 2 for x in data]
    mean_ssd_value = calculate_mean(ssd_values)
    ssd = float(sum(val for val in ssd_values if val <= 5 * mean_ssd_value))
    pvar = ssd / float(n - degrees_of_freedom)
    return pvar ** 0.5


def get_most_common_item(data):
    """
    Returns the most common item of a list. The list may contain tuples
    """
    assert isinstance(data, list)
    return max(set(data), key=data.count)


assert get_most_common_item([0]) == 0
assert get_most_common_item([0, 1, 1]) == 1
assert get_most_common_item([0, 1, 1, 0, 100, 0]) == 0
assert get_most_common_item([(0, 1), (0, 2), (0, 2), (0, 3)]) == (0, 2)
assert get_most_common_item([(0, 1), (0, 2), (0, 2), (0, 3), (0, 1), (0, 1)]) == (0, 1)
