========================================================================
    CONSOLE APPLICATION : Multi_Objective_Genetic_Algorithm Project Overview
========================================================================

NOTE: The result is an individual whose gender dominant attributes are successful, and whose
gender recessive attributes are also successful - they are good at gender stereotypical roles
-e.g. because the individuals father is good at hunting and the mother is good at foraging, the 
individual inherits the attributes from the parent.

This simple C++ program computes a multi object hierarchical genetic algorrithm. 
The heirarchy gene is simply male or female and this is implicit in the structure
of the program (2 lists). Each pair of parents produce 1 male and 1 female offspring.

The individuals are awarded greater breeding rights if they climb highest up the hill 
surface - but there is one catch - the females are walking on a different hill
surface to the males.

There are also genes shared between the males and females, for example these could be
genes for health, disease immunity, etc, they are also a seperate surface, but they
are summed into the fitness for each gender.


ADDED: There is now an experiement using gene coding ( GATC type coding for real genes).
This experiement is not in the stage of breeding, I have simply added some functionality
to support this ... the code will report back the proteins from a dna strand.

ADDED: some new objects for GA have been added, these will eventually be more generic 
and usable, however still in early stages.


