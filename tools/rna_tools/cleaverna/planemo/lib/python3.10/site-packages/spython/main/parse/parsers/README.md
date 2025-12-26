# Parsers

A parser class is intended to read in a container recipe file, and parse
sections into a spython.main.recipe Recipe object. To create a new subclass
of parser, you can copy one of the current (Docker or Singularity) as an
example, and keep in mind the following:

 - The base class, `ParserBase` in [base.py](base.py) has already added an instantiated (and empty) Recipe() for the subclass to interact with (fill with content).
 - The subclass is encouraged to define the name (self.name) attribute for printing to the user.
 - The subclass should take the input file as an argument to pass the the ParserBase, which will handle reading in lines to a list self.lines.
 - The subclass should have a main method, parse, that when called will read the input file and populate the recipe (and return it to the user).
