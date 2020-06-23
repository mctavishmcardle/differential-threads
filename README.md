# Differential Threads Reference

This repo contains:
* A list of the various possible differential thread combinations, among common
  UTS & ISO threads, along with the properties of those combinations:
  * The effective pitch & TPI
  * The clearance between the threads, if any, for combinations that could nest,
    rather than having to be formed on opposing ends of a single rod
* A script for generating that list

This is intended as a quick reference - a starting point for more rigorous
investigation of the possibilities of a given combination of threads. As it is
implemented using, for the most part, floating point arithmetic & leaves aside
entirely the question of tolerance, the values here cannot be considered
definitive for engineering purposes.
