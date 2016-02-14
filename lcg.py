import random as rnd
import numpy as np
from itertools import islice

"""
@Author: Tony Poerio
@email:  adp59@pitt.ed
University of Pittsburgh
Spring 2016
CS1538 - Simulation
Assignment #2 - Linear Congruential Random Number Generator

This python file will create a series of random number generators, and allow a user to specify
inputs for statistical tests, testing how truly random and indpendent each output is, exactly.

Tests performed:
    * Chi-squared for Uniformity
    * Kolmogorov-Smirnov Test for Uniformity
    * Runs Test for Independence
    * Autocorrelation test for independence at gaps:  2,3,5,50


*************************
WARNING:  This program will output test statistics for any input. But some of test suite is hard-coded with with
critical values for N=10_000, because that was the specified expected input for this project. These can be changed
in-situ below, if needed, any test result will be calculated as expected.
*************************

"""


#
# Formula to implement X_(i+1) = (aX + c) mod m
#       where a, c, and m are constants we choose


def main():
    # CONTROL FLOW
    print "UNIVERISTY OF PITTSBURGH - SPRING 2016:: CS1538, Assignment #2"
    print "--------------------------------------------------------------"

    # Get input, strip newlines.
    test_selection = ""
    while (test_selection != "q" ):
        select_test()
        test_selection = raw_input("Selection > ").strip()
        if test_selection == "q":
            exit()

        select_number_of_observations()
        number_observations = raw_input("Selection > ").strip()
        number_observations = int(number_observations)

        # If use selects python rand function,
        # create output file, and run the battery of tests
        if int(test_selection) == 1:
            python_rand( number_observations )
            run_test_suite(test_selection, number_observations)

        # If use selects LCG function,
        # create output file, and run the battery of tests
        elif int(test_selection) == 2:
            generate_lcg( number_observations )
            run_test_suite(test_selection, number_observations)



        # If user selects LCG with RANDU settings
        # create output file, and run the battery of tests
        elif int(test_selection) == 3:
            generate_lcg_RANDU( number_observations )
            run_test_suite(test_selection, number_observations)
        else:
            print "Please select a number from 1 to 3."



    #   - select a function to generate random numbers
    #   - how many random numbers to generate
          # >>> HOW MANY?

    # get num observations, set to 100 for now


    # Generate output files
    # python_rand( number_observations )
    # generate_lcg( number_observations )
    # generate_lcg_RANDU( number_observations )

    """
    # divide our output values in 10 equal subdivisions and run chi-square test
    data_points = divide_RNG_data_into_10_equal_subdivisions_and_count("lgc_output.txt")
    result = chi_square_uniformity_test(data_points, 0, number_observations)
    print "Result of chi-square = " + str(result)

    # get first 100 values from sample and run kolmogorov-smirnov test
    first_100_values = collect_first_100_samples_in_data_set("lgc_output.txt")
    first_100_values.sort()
    ks_result = kolmogorov_smirnov_test(first_100_values,1,100)

    # perform a runs test
    runs_test_result = runs_test_for_independence("lgc_output.txt", number_observations )
    print "Runs Test Result Z-Score: " + str(runs_test_result)

    # perform an autocorrelation test
    auto_test_result = autocorrelation_tests("lgc_output.txt", number_observations, 10 )

    """

    #   - select a statistical test to run
    #   - save a stream of generate numbers to a file, rather than let them scroll across the screen

    # THREE GENERAL FUNCTION SETTINGS
    # Same seed for all cases (seed=123456789)
    # x 1.  Standard random number generator in Python
    # x 2.  LCG Implementation
    #         o Where:  a=101427; c=321, m=(2**16)
    #         o Obtain each number in U[0,1) by diving X_i by m
    # x 3.  LCG with RANDU
    #         o Where: a=65539; c=0; m=(2**31)
    #         o Again, obtain each number in U[0,1) by diving X_i by m



#########################
###   RNG FUNCTIONS   ###
#########################

def python_rand( num_iterations ):
     """
     Run the built-in python random number generator and output a number of data points
     specified by the user to a file
     :param num_iterations:  The number of data points to write to file
     :return: void
     """
     # Initialize seed value
     x_value = 123456789.0    # Our seed, or X_0 = 123456789
     rnd.seed(x_value)

     # counter for how many iterations we've run
     counter = 0

     # Open a file for output
     outFile = open("py_random_output.txt", "wb")

     # Perform number of iterations requested by user
     while counter < num_iterations:
        x_value = rnd.random()
        # Write to file
        writeValue = str(x_value)
        outFile.write(writeValue + "\n")
        counter = counter + 1

     outFile.close()
     print("Successfully stored %d random numbers in file named: 'py_random_output.txt'.", num_iterations)


def generate_lcg( num_iterations ):
    """
    LCG - generates as many random numbers as requested by user, using a Linear Congruential Generator
    LCG uses the formula: X_(i+1) = (aX_i + c) mod m
    :param num_iterations: int - the number of random numbers requested
    :return: void
    """
    # Initialize variables
    x_value = 123456789.0    # Our seed, or X_0 = 123456789
    a = 101427               # Our "a" base value
    c = 321                  # Our "c" base value
    m = (2 ** 16)            # Our "m" base value

    # counter for how many iterations we've run
    counter = 0

    # Open a file for output
    outFile = open("lgc_output.txt", "wb")

    #Perfom number of iterations requested by user
    while counter < num_iterations:
        # Store value of each iteration
        x_value = (a * x_value + c) % m

        #Obtain each number in U[0,1) by diving X_i by m
        writeValue = str(x_value/m)

        # write to output file
        outFile.write(writeValue + "\n")
        # print "num: " + " " + str(counter) +":: " + str(x_value)

        counter = counter+1

    outFile.close()
    print("Successfully stored " + str(num_iterations) + " random numbers in file named: 'lgc_output.txt'.")


def generate_lcg_RANDU( num_iterations ):
    """
    LCG RANDU- generates as many random numbers as requested by user, using a Linear Congruential Generator
    LCG uses the formula: X_(i+1) = (aX_i + c) mod m.

    This LCG uses the RANDU initial setting, a=65539; c=0; m=2^31.
    RANDU is known to have an issue: its values fall into 15 parallel 2D planes.
    So while its pseudo-randomness is enough for some applications. It's not great.
    Not crypto strength by any means.

    :param num_iterations: int - the number of random numbers requested
    :return: void
    """
    # Initialize variables
    x_value = 123456789.0    # Our seed, or X_0 = 123456789
    a = 65539                # Our "a" base value
    c = 0                    # Our "c" base value
    m = (2 ** 31)            # Our "m" base value

    # counter for how many iterations we've run
    counter = 0

    # Open a file for output
    outFile = open("lgc_RANDU_output.txt", "wb")

    #Perfom number of iterations requested by user
    while counter < num_iterations:
        # Store value of each iteration
        x_value = (a * x_value + c) % m

        #Obtain each number in U[0,1) by diving X_i by m
        writeValue = str(x_value/m)

        # write to output file
        outFile.write(writeValue + "\n")
        # print "num: " + " " + str(counter) +":: " + str(x_value)

        counter = counter+1

    outFile.close()
    print "Successfully stored " + str(num_iterations) + " random numbers in file named: 'lgc_RANDU_output.txt'."


######################
#### STATS TESTS #####
######################
    # STATISTICAL TESTS
    # Check for uniformity at 80%, 90%, and 95% level. Note that some tests are one-sided, others two sided
    # x 1. Chi-Square Frequency Test for Uniformity
    #      - Collect 10,000 numbers per generation method
    #      - Sub-divide[0.1) into 10 equal subdivisions
    # x 2. Kolmogorov-Smirnov Test for uniformity
    #      - Since K-S Test works better with a smaller set of numbers, you may use the first 100
    #        out fo the 10,000 that you generated for the Chi-Square Frequency Test
    # x 3. Run Test for Independence
    #      - Use all 10,000 numbers per method for this test
    # x 4. Autocorrelations test for independence
    #      - Use all 10,000 numbers per method
    #      - Set of tests, but consider different parameterizations.
    #      - Experiment with diff gap sizes (param 1 in our notes): 2,3,5,50.
    #      - But don't consider different starting points (parameter i in our notes)
    #      - Simply choose one and specify what was used, in the report
    #      - Check these links to help implement: http://www.itl.nist.gov/div898/handbook/eda/section3/eda35c.htm
    #      - http://www.aritzhaupt.com/resource/autocorrelation/


def chi_square_uniformity_test( data_set, confidence_level, num_samples ):
    """
    Null hypothesis:  Our numbers distributed uniformly on the interval [0, 1).

    This function uses the chi-square test for uniformity to determine whether our numbers
    are uniformly distributed on the interval [0,1).

    Formula is: "sum[ (observed-val - expected-val)^2 / expected val ], from 0 to num_samples"
    This gives us a number which we can test against a chi-square value table.

    Also need to know, degrees of freedom:  df=num_samples-1
    :param data_set: the data_set, must be a dictionary with 10 intervals.
                     Use return value from  @divide_RNG_data_into_10_equal_subdivisions_and_count
    :param confidence_level: confidence level we are testing at
    :param num_samples: number of data points
    :return: A chi-squared value
    """
    # This is our test statistic, this will be an accumulated value, as we loop through the data set
    chi_sq_value = 0.0
    degrees_of_freedom = num_samples - 1

    # We're doing 10 equal subdivisions, so need to divide our number samples by 10,
    # Assuming uniform distribution, to get an expected value. All values should be same
    # If our distro is actually uniform.
    expected_val = num_samples/10.0


    # Loop through a dictionary and get every count
    # The observed value is going to be our count at each key, and then we can do chi-square
    for observed_val in data_set:
        # print "Observed value is: " + observed_val
        chi_sq_value += ( pow((expected_val - data_set[observed_val]), 2)/expected_val )

    # Coming out of this loop, we'll have a chi-squared test statistic
    # Now we just need to do a lookup to see if it's valid
    return chi_sq_value



def kolmogorov_smirnov_test( data_set, confidence_level, num_samples ):
    """
    Kolmogorov-Smirnov test for uniform distribution of Random numbers
    :param data_set: The set of data to analyze. Should be floating point numbers [0,1) in a .txt file
    :param confidence_level: with how much confidence should we test?
    :param num_samples: number of samples to analyze
    :return: test statistic
    """
    # Step 1:  Rank data from smallest to largest, such that:
    # R(1) <= R(2) <= R(3) ... <= R(i)
    data_set.sort()

    # Step 2: Computer D+ and D-
    # D+ = max(i/N - R(i))
    d_plus = get_d_plus_value_for_KS_TEST(data_set, num_samples)
    print "D+ VALUE ="+str(d_plus)

    # D- = max(R(i) - (i -1)/n)
    d_minus = get_d_minus_value_for_KS_TEST(data_set, num_samples)
    print "D- VALUE="+str(d_minus)

    # Step 3:  Computer D = max(D+,D-)
    d_value = max(d_plus, d_minus)
    print "D VALUE (max): "+str(d_value)

    # Step 4: Determine critical value, using table
    # Step 5: Accept or reject Null hypothesis
    return d_value


def runs_test_for_independence( data_file, num_samples ):
    """
    Perform a runs test for independence for a set of random numbers
    :param data_set: A data set with 10,000 samples - should be a list, with all numbers floats
    :param num_samples: The number of samples to test
    :return: z-test statistic for our data set
    """
    # Note, every sequence begins and ends with "NO EVENT"
    # Run = Succession of similar events, followed by a different event
    # Run Length = Number of events that occur in the run
    # Number of runs = number of "runs" total
    # Two concerns: Num runs, length of runs
    # We're looking for:  Runs of larger and smaller numbers (increasing or decreasing)
    runLengths = { }       # The size of each run, in order
    numRuns = 0            # The number of runs overall
    runDirection = "none"  # We'll use "none", "up", and "down" to keep track

    with open( data_file, "r" ) as f:
        data_points = f.readlines()


   # for value in data_points:
    for value in range(0, (len(data_points)-1) ):

        # DEBUG
        thisValue = float(data_points[value])
        nextValue = float(data_points[value+1])
        # print "data_points[value] ::"+ str(thisValue)
        # print "dat_points[value+1] ::"+ str(nextValue)

        # If no change in direction we'll ignore
        if thisValue == nextValue:
            numRuns = numRuns             # numRums overall, doesn't change
            runDirection = runDirection   # runDirection doesn't change

        # Check if we have a NEW run, going UP
        elif thisValue < nextValue and runDirection != "up":
            numRuns = numRuns + 1         # We have a NEW run
            runDirection = "up"           # We have a NEW run direction
            runLengths[numRuns] = 1       # We have a NEW key in our dictionary, with value=1

        # Check if we have a CONTINUING run, going UP
        elif thisValue < nextValue and runDirection == "up":
            runLengths[numRuns] += 1      # increment the run length in our dictionary for current run
                                          # NumRuns doesn't change
                                          # runDirection doesn't change

        # Check if we have NEW run, going DOWN
        elif thisValue > nextValue and runDirection != "down":
            numRuns = numRuns + 1          # We have a NEW run
            runDirection = "down"          # We have a NEW run direction
            runLengths[numRuns] = 1        # We have a NEW key in our dictionary, with value=1

        # Check if we have a CONTINUING run, going DOWN
        elif thisValue > nextValue and runDirection == "down":
            runLengths[numRuns] += 1      # increment the run length in our dictionary for current run
                                          # NumRuns doesn't change
                                          # runDirection doesn't change

    # Leaving this loop, we should have a dictionary with our run numbers mapped to their lengths
    # We should also have a the number of runs

    # Now, calculate mean:  Mean = (2N-1)/3
    mean =  ( (2*num_samples - 1) / 3 )

    # And variance:  Variance = (16(N) - 29) / 90
    variance = ( ( 16*num_samples - 29) / 90 )

    # And we can use the mean & variance to calculate the Z-Test statistic
    z_statistic = ( (numRuns - mean) / np.sqrt(variance) )

    print "Number of runs: " + str(numRuns)

    return z_statistic


def autocorrelation_tests( data_file, num_samples, gap_sequence ):
    """
    This function performs an autocorrelation test
    :param data_file:  The data set to analyze. Should be 10,000 samples, all floats, in a .txt file
    :param num_samples:  The number of samples to analyize
    :return: The test result
    """

    # Get necesssary values, store in "data_points" list
    with open( data_file, "r" ) as f:
         data_points = f.readlines()

    # Sort data set
    # data_points.sort()

    little_m = gap_sequence    # The space between the numbers being tested
    ###### CHANGE BACK TO 0 #########
    start_index = 0            # The number in the sequence we start with
    ########### end #################
    big_n = num_samples        # the number of numbers generated in a sequence
    big_m = 0.0                  # Largest number such that i + (M+1)m <= N

    # Determine correct M value
    while (big_m + 1) < ( (big_n - start_index)/little_m ) :
        big_m = big_m + 1

    # print "Final value for big_m: " + str(big_m)

    one_over_m_plus_one = ( 1.0/(big_m + 1.0 ) )
    rho_hat = 0.0
    sum_of_rho_hat = 0.0


    # Get every m'th element in the data_set
    every_m_element = data_points[0::gap_sequence]

    # print "Length: " + str(len(every_m_element))
    # Get the sum of rho_hat
    for value in range(0, (len(every_m_element)-1) ):
        thisValue = float(every_m_element[value])
        nextValue = float(every_m_element[value+1])
        # print "Autocorrelation: Ki   :" + str(thisValue)
        # print "Autocorrelation: Ki+1 :" + str(nextValue)
        sum_of_rho_hat = sum_of_rho_hat + (thisValue * nextValue)
        # print "Sum of rho hat: " + str(sum_of_rho_hat)

    # Subtract 0.25
    sum_of_rho_hat = (one_over_m_plus_one * sum_of_rho_hat) - 0.25

    variance_of_rho =  np.sqrt( (13*big_m + 7 )) / (12*(big_m + 1))

    z_statistic = sum_of_rho_hat / variance_of_rho

    # print "Z-Score for autocorrelation is: " + str(z_statistic)
    return z_statistic



##############################
##### Significance Tests #####
##############################

def chi_sq_significance_test( chi_sq, signif_level):
    """
    Performs a significance test for df=10000, based on values calculated at:
    https://www.swogstat.org/stat/public/chisq_calculator.htm
    :param chi_sq:  Chi-sq value to test
    :param signif_level: Level of significance we are testing: 0.80, 0.90, or 0.95
    :return: message stating whether we accept or reject null
    """
    result = "FAIL TO REJECT null hypothesis"
    crit_value = 0.0
    if signif_level == 0.8:
        crit_value = 10118.8246
    elif signif_level == 0.90:
        crit_value = 10181.6616
    elif signif_level == 0.95:
        crit_value = 10233.7489
    else:
        print "**Invalid Significance Level for Chi Sq***"

    if chi_sq > crit_value:
        result = "REJECT null hypothesis"

    print "Print Significance Level: " + str(signif_level)
    print "Chi Sq: " + str(chi_sq)
    print "Crit Value: " + str(crit_value)
    print "Result is: " + result
    print "...................................."

    return result

def ks_significance_test( d_statistic, num_observations, alpha_level ):
    """
    Perform Significance test for Kolmogorov-Smirnov
    Uses formulas from table A.7:  Discrete-Event System Simulation, by Banks and Carson, 1984
    :param d_statistic: The d-value we are testing
    :param num_observations: The number of observations in our data set
    :param alpha_level: The level of significance we are testing
    :return: result -- accept or reject
    """
    result = "FAIL TO REJECT null hypothesis"
    critical_value = 0


    if alpha_level == 0.1:
        critical_value = 1.22/np.sqrt(num_observations)
    elif alpha_level == 0.05:
        critical_value = 1.36/np.sqrt(num_observations)
    elif alpha_level == 0.01:
        critical_value = 1.63/np.sqrt(num_observations)
    else:
        print "Invalid alpha level for KS test. Must be: 0.1, 0.05, or 0.01"

    if d_statistic > critical_value:
        result = "REJECT null hypothesis"
    print "Alpha Level is: " + str(alpha_level)
    print "D_statistic is: " + str(d_statistic)
    print "Critical value is: " + str(critical_value)
    print "Result is: " + result
    print "............................"

    return result

def z_score_lookup( z_score, significance_level, two_sided=True):
    """
    Performs a two-sided z-score lookup, for 0.8, 0.9, or 0.95 level of significance
    :param z_score: Z score to test
    :param significance_level: Significance level
    :return: String detailing our result
    """

    result = "FAIL TO REJECT null hypothesis"
    critical_value = 0.0
    confidence_80 = 1.282
    confidence_90 = 1.645
    confidence_95 = 1.96
    confidence_99 = 2.576

    # Assign confidence interval z-scores to our crit value
    if significance_level == 0.8:
        critical_value = confidence_80
    elif significance_level == 0.9:
        critical_value = confidence_90
    elif significance_level == 0.95:
        critical_value = confidence_95
    else:
        print "Invalid significance level for z-lookup. Must be: 0.8, 0.9, or 0.95"

    # Need to adjust intervals if the test is one sided
    if two_sided == False:
        if critical_value == confidence_80:
            critical_value = 0.8416
        elif critical_value == confidence_90:
            critical_value = 1.282
        elif critical_value == confidence_95:
            critical_value = 1.645



    neg_crit_value = critical_value * (-1.0)

    #if z_score < 0:
     #  z_score = z_score * (-1)

    if ( two_sided and ( z_score <= neg_crit_value) or (critical_value <= z_score ) ):
        result = "REJECT null hypothesis"

    if ( not two_sided and z_score >= critical_value or z_score <= neg_crit_value):
        result = "REJECT null hypothesis"

    print "Z score is: " + str(z_score)
    print "Significance level is: " + str(significance_level)
    print "Critical value is: " +str(critical_value)
    print "Running two sided z-score lookup? -->" + str(two_sided)
    print ""
    print "Result is: " + result
    print "....................................."

    return result



#######################
### Helper Methods ####
#######################

def collect_first_100_samples_in_data_set( data_file ):
    """
    Takes a data file, with real number data points between [0,1) reads the first 100 values,
    then adds them to a dictionary as our return value
    :param data_file: A string - the name of the file to read in our current directory
    :return: A dictionary containing the first 100 values as floats
    """

    first_100_vals_as_FLOATS = []
    # grabs first 100 files, as strings with newline endpoints
    with open( data_file, "r" ) as f:
        first_100_vals_as_STRINGS = list(islice(f, 100))

    # transform all values to floats
    for val in first_100_vals_as_STRINGS:
        val = float(val)
        first_100_vals_as_FLOATS.append(val)

    return first_100_vals_as_FLOATS


def divide_RNG_data_into_10_equal_subdivisions_and_count( data_file ):
    """
    Takes a path to a data file in the current directory.
    Returns a dictionary with keys 1-10, values=num instances in each of
    10 equal intervals from range: [0, 1).
    The function counts how many data points are in each interval, and gives us
    a dictionary so we can manipulate this data more easily, based on count by index.

    :param data_file: Must be in current directory. Pass in the string name.
    :return: A dictionary with counts of how many occurrences our data had for each
    of 10 equal intervals between [0, 1). (Divided into 10ths)
    """
    # For each of our uniformity tests, need to divide our data points in 10 equal subdivisions
    subdivisions = {  "1":  0,
                      "2":  0,
                      "3":  0,
                      "4":  0,
                      "5":  0,
                      "6":  0,
                      "7":  0,
                      "8":  0,
                      "9":  0,
                      "10": 0   }
    with open(data_file, "r") as f:
        # data points is a list containing all numbers we've read in.
        data_points = f.readlines()

    # Loop through our data points and count number of data points in each subdivision
    # Divide by tenths, from 0.0 to 1.0.
    for num in data_points:
        num = float(num)
        if num < 0.1:
            subdivisions["1"] += 1
        elif num < 0.2:
            subdivisions["2"] += 1
        elif num < 0.3:
            subdivisions["3"] += 1
        elif num < 0.4:
            subdivisions["4"] += 1
        elif num < 0.5:
            subdivisions["5"] += 1
        elif num < 0.6:
            subdivisions["6"] += 1
        elif num < 0.7:
            subdivisions["7"] += 1
        elif num < 0.8:
            subdivisions["8"] += 1
        elif num < 0.9:
            subdivisions["9"] += 1
        elif num < 1.0:
            subdivisions["10"] += 1

    return subdivisions


def get_d_plus_value_for_KS_TEST( data_set, num_samples ):
    """
    Finds the D+ value for a KS test
    :param data_set: 100 values, must be a list of floats
    :return: the D-+Statistic for our data set
    """
    # D+ = max(i/N - R(i))
    d_plus_max = 0
    value_rank_i = 1

    # iterate through data set
    for value in data_set:
        # Do each D+ calculation, store it
        d_plus_i_value = ( (value_rank_i/num_samples) - value )

        # Check if it is highest D+ value yet
        if d_plus_i_value > d_plus_max:
            d_plus_max = d_plus_i_value

        # increment our "i" value
        value_rank_i = value_rank_i + 1

    # coming out of this loop, D+ = highest D+ value
    return d_plus_max


def get_d_minus_value_for_KS_TEST( data_set, num_samples ):
    """
    Finds the D- value for a KS test
    :param data_set: 100 values, must be a list of floats
    :return: the D- Statistic for our data set
    """
    # D- = max(R(i) - (i -1)/n)
    d_minus_max = 0
    value_rank_i = 1.0

    # iterate through data set
    for value in data_set:
        # Do each D+ calculation, store it
        substraction_value = ( (value_rank_i - 1.0)/num_samples )
        d_minus_i_value = value - substraction_value

        # Check if it is highest D+ value yet
        if d_minus_i_value > d_minus_max:
            d_minus_max = d_minus_i_value

        # increment our "i" value
        value_rank_i = value_rank_i + 1

    # coming out of this loop, D+ = highest D+ value
    return d_minus_max

def select_test():
    """
    Command line prompt for selecting a test
    :return: void - prints a prompt to command line
    """
    print "Please select a method for generating random numbers: "
    print " 1. Python's Random Function"
    print " 2. Linear Congruential Generator "
    print " 3. Linear Congruential Generator with RANDU initial settings"
    print ""
    print "      (or type 'q' to quit)"
    print ""

def select_number_of_observations():
    """
    Command line prompt to select the number of observations for a given test
    :return: void - prints a prompt to command line
    """
    print "How many observations should we perform?"


def run_test_suite( test_selection, number_observations ):
    """
    Runs all of our test suites and prints output to the screen
    :param test_selection:  an int - 1,2, or 3. Corresponds to test selected.
    :param number_observations: the number of data points to test
    :return: void - prints to command line
    """
    input_file = ""
    test_name = ""
    test_selection = int(test_selection)
    if test_selection == 1:
        input_file = "py_random_output.txt"
        test_name = "PYTHON BUILT-IN RAND"

    elif test_selection == 2:
        input_file = "lgc_output.txt"
        test_name = "LINEAR CONGRUENTIAL GENERATOR"

    elif test_selection == 3:
        input_file = "lgc_RANDU_output.txt"
        test_name = "LGC with RANDU initial settings"
    else:
        print "Invalid input. Please try again."

    print ""
    print ""
    print "TEST SUITE FOR:  %s " % (test_name)
    print "======================================"

    # divide our output values in 10 equal subdivisions and run chi-square test
    print "---------CHI-SQ_TEST-----------"
    data_points = divide_RNG_data_into_10_equal_subdivisions_and_count(input_file)
    chi_sq_result = chi_square_uniformity_test(data_points, 0, number_observations)
    chi_sq_significance_test( chi_sq_result, 0.8 )
    chi_sq_significance_test( chi_sq_result, 0.9 )
    chi_sq_significance_test( chi_sq_result, 0.95 )

    print ""

    # get first 100 values from sample and run kolmogorov-smirnov test
    print "---------KS_TEST-----------"
    first_100_values = collect_first_100_samples_in_data_set(input_file)
    first_100_values.sort()
    ks_result = kolmogorov_smirnov_test(first_100_values,1,100)
    ks_significance_test(ks_result,100, 0.1)
    ks_significance_test(ks_result,100, 0.05)
    ks_significance_test(ks_result,100, 0.01)
    print "Kolmogorov-Smirnov Test Result for D-Value: " + str(ks_result)
    print ""

    # perform a runs test
    print "---------RUNS_TEST-----------"
    runs_test_result = runs_test_for_independence( input_file, number_observations )
    print "Runs Test Result Z-Score: " + str(runs_test_result)
    print ""
    z_score_lookup(runs_test_result, 0.8, two_sided=True)
    z_score_lookup(runs_test_result, 0.9, two_sided=True)
    z_score_lookup(runs_test_result, 0.95, two_sided=True)
    print ""

    # perform an autocorrelation test
    print "---------AUTO-CORRELATION_TESTS-----------"
    auto_test_result = autocorrelation_tests(input_file, number_observations, 2 )
    print "==== Auto-correlation Test Result for GAP SIZE=2:  " + str(auto_test_result)

    z_score_lookup(auto_test_result, 0.8, two_sided=True)
    z_score_lookup(auto_test_result, 0.9, two_sided=True)
    z_score_lookup(auto_test_result, 0.95, two_sided=True)
    print ""
    print "       ===== END GAPSIZE=2 ====="
    print ""
    print ""
    auto_test_result = autocorrelation_tests(input_file, number_observations, 3 )
    print "=== Auto-correlation Test Result for GAP SIZE=3: " + str(auto_test_result)
    z_score_lookup(auto_test_result, 0.8, two_sided=True)
    z_score_lookup(auto_test_result, 0.9, two_sided=True)
    z_score_lookup(auto_test_result, 0.95, two_sided=True)
    print ""
    print "       ===== END GAPSIZE=3 ====="
    print ""
    print ""
    auto_test_result = autocorrelation_tests(input_file, number_observations, 5 )
    print " === Auto-correlation Test Result for GAP SIZE=5: " + str(auto_test_result)
    z_score_lookup(auto_test_result, 0.8, two_sided=True)
    z_score_lookup(auto_test_result, 0.9, two_sided=True)
    z_score_lookup(auto_test_result, 0.95, two_sided=True)
    print ""
    print "       ===== END GAPSIZE=5 ====="
    print ""
    print ""
    auto_test_result = autocorrelation_tests(input_file, number_observations, 50 )
    print "Auto-correlation Test Result for GAP SIZE=50: " + str(auto_test_result)
    z_score_lookup(auto_test_result, 0.8, two_sided=True)
    z_score_lookup(auto_test_result, 0.9, two_sided=True)
    z_score_lookup(auto_test_result, 0.95, two_sided=True)
    print ""
    print "       ===== END GAPSIZE=50 ====="
    print ""
    print ""
    print ""
    print "=-=-=-======= END TEST SUITE ======-=-=-=-="
    print ""

if __name__ == "__main__":
    main()




""" Write-Up / Analysis

Once you have done all of your tests, you must look over and analyze your results and present them in a well-written,
well-formatted report. Your report should include a discussion of the following:

What is the random library function that your compiler supplied for you?
Summarize the outcomes of the statistical tests for each RNG method in a formatted, easy to understand table.
What test(s) did the method "pass"? Be specific about the condition (i.e., at what level of significance).
For each method, prior to generating the numbers and running the statistical tests,
 do you expect it to work well or poorly? Explain.
Looking over some of the generated numbers for each method (but before running the statistical tests),
do you think they look sufficiently random? Did the outcome of the statistical test surprise you?
Discuss whether you think the set of experiments you did for this assignment is sufficient.
If so, argue why that is the case. If not, explain what additional test(s) or modification(s) to the methodology you'd perform.
Submission and Grading:

The assignment is due Wednesday, February 17 by 11:59 pm.

Zip your source code files, your write-up, and your README into one zip file and upload it to CourseWeb, in the Homework 2 location.
For more advice on submitting your assignment, see the Assignments section of the Tips for Success page.
"""""


