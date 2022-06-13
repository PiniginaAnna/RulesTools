# RulesTools

Tools for extracting different types of reaction rules and filtering wrong reactions


## Installing

```
pip install git+https://gitlab.cimm.site/pinigina/rulestools
```


## Usage

Extract individual reaction rule

```python
from RulesTools.transformations import CreateRule

rules_creator = CreateRule(rules_from_multistage_reaction=False, environment_atoms_number=1, 
                           rule_with_functional_groups=True, functional_groups_list=groups_list, include_rings=True,
                           keep_reagents=False, keep_meta=False, as_query=True, keep_atom_info='reaction_center')
rules = rules_creator(reaction)
```

Extract reaction rules from the reaction database

```python
from RulesTools.database_processing import reaction_database_processing
from RulesTools.transformations import *
from RulesTools.filters import *

transformations = [DeleteSmallMolecules(number_of_atoms=4), ReverseReaction(),
                   CreateRule(rules_from_multistage_reaction=False, environment_atoms_number=1, 
                              rule_with_functional_groups=True, functional_groups_list=groups_list, include_rings=True,
                              keep_reagents=False, keep_meta=False, as_query=True, keep_atom_info='reaction_center')]

reaction_database_processing('path_to_reaction_database.rdf', transformations=transformations, save_only_unique=True)
```
