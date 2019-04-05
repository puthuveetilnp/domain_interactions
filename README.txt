To calculate a domain pair frequency list, first you need:

1. A paired list of bait and prey proteins. See E.coli_paired.csv for example
2. A list from Uniprot with the domains (PFAM or InterPro) of the bait proteins. See E.coli_bait.csv for example
3. A list from Uniprot with the domains (PFAM or InterPro) of the bait proteins. See E.coli_prey.csv

Next, open and run the file "find_domain_pairs.py"

The program will first prompt for the name of the organism you are interested in. If you're interested
in E.coli domains, then type E.coli.

Next, the program will prompt for the name of the pair list file, then bait list file, and the prey list.

Lastly, the program will ask if you are using PFAM ids or InterPro ids.
