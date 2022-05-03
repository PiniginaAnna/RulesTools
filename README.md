# RulesTools

Tools for extracting different types of reaction rules


## Installing

```
Give the example
```


## Usage

```python

transformations = [DeleteSmallMolecules(number_of_atoms=4), ReverseReaction(),
                   CreateRule(rules_from_multistage_reaction=False, environment_atoms_number=1, include_rings=False,
                                rule_with_functional_groups=False, functional_groups_list=groups_list, as_query=True,
                                keep_atom_info='reaction_center', keep_reagents=False, keep_meta=False)]
result_reactions = apply_transformations(transformations, reaction)
result_reactions
```

The results contain the mapped reactions and confidence scores:

```python

reaction_rules_extraction(reaction_database_file, transformations, directory_name)

```
