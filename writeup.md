**Anthony Poerio - (adp59@pitt.edu)**
**CS1538 - Assignment #2, Writeup**
**February 13, 2016**

## Random Number Generation
Random numbers… because deterministic… goal is to get a long period between repeats… and have no patterns — to be uniformly distributed.

Null hypothesis… 

Linear Congruential Generators are …. which… 

In this assignment, three (3) total Linear Congruential Generators:
	1. Python’s Built in Random Function
	2. A Linear Congruential Generator (LCG) with Seed value: 
	3. A Linear Congruential Generator with the RANDU initial settings

In this report, we will run a suite of tests on the output of each generator, attempting to determine just how random each function really is. All random functions except the RANDU are seeded with the value: 123456789. To replicate each test described here, please run and see the source code in: lcgy.py, which is attached. 

## Python Random 
The random library function supplied by the Python programming language is called by invoking random(), which calls Pseudo-Random Number generation algorithm known as the  “Mersenne Twister” (PyDocs). The Mersenne Twister makes use of very large prime numbers to  is known to have a period of **(2^(19937)-1)**. It is known to pass many statistical randomness tests (Mersenne Twister Wiki), And the results of my testing we as follows: 

Large Chi square… value means there is strong relationship between our variables.


We can reject independence with 80% certainty, and 90% certainty… but not 95% certainty

Runs and Autocorrelation test Independence
https://www.swogstat.org/stat/public/chisq_calculator.htm

| Test Name                            	| Sample Size      	| Significance Level 	| Critical Value 	| Test Statistic Found 	|             Result             	|
|--------------------------------------	|------------------	|--------------------	|----------------	|---------------------:	|:------------------------------:	|
| Chi-Square                           	| 10000            	| 0.80               	| 118.5          	| 9.754                	| FAIL TO REJECT Null Hypothesis 	|
| Chi-Square                           	| 10000            	| 0.90               	|                	| 9.754                	| FAIL TO REJECT Null Hypothesis 	|
| Chi-Square                           	| 10000            	| 0.95               	|                	| 9.754                	| FAIL TO REJECT Null Hypothesis 	|
| Kolmogorov-Smirnov                   	| first 100 values 	|                    	|                	|                      	|                                	|
| Runs Test                            	| 10000            	|                    	|                	|                      	|                                	|
| Autocorrelation Test with GapSize=2  	| 10000            	|                    	|                	|                      	|                                	|
| Autocorrelation Test with GapSize=3  	| 10000            	|                    	|                	|                      	|                                	|
| Autocorrelation Test with GapSize=5  	| 10000            	|                    	|                	|                      	|                                	|
| Autocorrelation Test with GapSize=50 	| 10000            	|                    	|                	|                      	|                                	|

## 

Once you have done all of your tests, you must look over and analyze your results and present them in a well-written, well-formatted report. Your report should include a discussion of the following:
	•	What is the random library function that your compiler supplied for you?
	•	Summarize the outcomes of the statistical tests for each RNG method in a formatted, easy to understand table. What test(s) did the method "pass"? Be specific about the condition (i.e., at what level of significance).
	•	For each method, prior to generating the numbers and running the statistical tests, do you expect it to work well or poorly? Explain.
	•	Looking over some of the generated numbers for each method (but before running the statistical tests), do you think they look sufficiently random? Did the outcome of the statistical test surprise you?
	•	Discuss whether you think the set of experiments you did for this assignment is sufficient. If so, argue why that is the case. If not, explain what additional test(s) or modification(s) to the methodology you'd perform.

## Works Cited
Python Random Library Docs:  https://docs.python.org/2/library/random.html 

Marsenne Twister:  https://en.wikipedia.org/wiki/Mersenne_Twister#cite_note-32 