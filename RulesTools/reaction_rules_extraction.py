from CGRtools import RDFRead, RDFWrite, SDFRead
from classes import *
from tqdm import tqdm
from collections import defaultdict
from pickle import dump
from os import makedirs


def apply_transformations(transformations, reaction):
    for transform in transformations:
        reaction = transform(reaction)
    return reaction


def reaction_rules_extraction(reaction_database_file_name, transformations, directory_name = './',
                              error_reactions_file_name = 'error_reactions', rules_file_name = 'rules'):

    unique_rules = defaultdict(lambda: [0, []])

    with RDFRead(reaction_database_file_name, indexable=True) as reactions:
        reactions.reset_index()

        for reaction in tqdm(reactions, total=len(reactions)):
            try:
                rules = apply_transformations(transformations, reaction)
            except Exception:
                with RDFWrite(f'{directory_name}/{error_reactions_file_name}.rdf', append=True) as errors_file:
                    errors_file.write(reaction)
            else:
                for rule in rules:
                    unique_rules[rule][0] += 1
                    unique_rules[rule][1].append(reaction.meta['Reaction_ID'])

    with RDFWrite(f'{directory_name}/{rules_file_name}.rdf') as rules_file:
        with open(f'{directory_name}/{rules_file_name}.pickle', 'ab') as rules_file2:
            for rule, info in unique_rules.items():
                rule.meta['number_of_reactions'] = info[0]
                rule.meta['reactions_id'] = ' '.join(info[1])
                rules_file.write(rule)
                dump(rule, rules_file2)


if __name__ == '__main__':

    with SDFRead('./data/groups.sdf', indexable=True) as groups:
        molecule_groups = list(groups)
    query_groups = []
    for group in molecule_groups:
        if isinstance(group, QueryContainer):
            query_groups.append(group)
        else:
            group.canonicalize()
            query_group = group.substructure(group, as_query=True)
            for atom_number in query_group.atoms_numbers:
                query_group.atom(atom_number).neighbors = None
                query_group.atom(atom_number).hybridization = None
                query_group.atom(atom_number).implicit_hydrogens = None
                query_group.atom(atom_number).ring_sizes = None
            query_groups.append(query_group)

    all_transformations = {
        'env_1':
            [CreateRule(rules_from_multistage_reaction=False, environment_atoms_number=1, include_rings=False,
                        rule_with_functional_groups=False, functional_groups_list=query_groups, as_query=True,
                        keep_atom_info='reaction_center', keep_reagents=False, keep_meta=False)],
        'env_1_rings':
            [CreateRule(rules_from_multistage_reaction=False, environment_atoms_number=1, include_rings=True,
                        rule_with_functional_groups=False, functional_groups_list=query_groups, as_query=True,
                        keep_atom_info='reaction_center', keep_reagents=False, keep_meta=False)],
        'env_1_groups':
            [CreateRule(rules_from_multistage_reaction=False, environment_atoms_number=1, include_rings=False,
                        rule_with_functional_groups=True, functional_groups_list=query_groups, as_query=True,
                        keep_atom_info='reaction_center', keep_reagents=False, keep_meta=False)],
        'env_1_groups_rings':
            [CreateRule(rules_from_multistage_reaction=False, environment_atoms_number=1, include_rings=True,
                        rule_with_functional_groups=True, functional_groups_list=query_groups, as_query=True,
                        keep_atom_info='reaction_center', keep_reagents=False, keep_meta=False)],
        'env_2':
            [CreateRule(rules_from_multistage_reaction=False, environment_atoms_number=2, include_rings=False,
                        rule_with_functional_groups=False, functional_groups_list=query_groups, as_query=True,
                        keep_atom_info='reaction_center', keep_reagents=False, keep_meta=False)],
        'groups_rings':
            [CreateRule(rules_from_multistage_reaction=False, environment_atoms_number=0, include_rings=True,
                        rule_with_functional_groups=True, functional_groups_list=query_groups, as_query=True,
                        keep_atom_info='reaction_center', keep_reagents=False, keep_meta=False)],
    }

    reaction_database_file = './data/short_USPTO_1000.rdf'
    all_results_directory = 'short_USPTO_1000_test'

    for transform_name, transformations in all_transformations.items():
        print(transform_name)
        directory_name = f'./{all_results_directory}/{transform_name}'
        makedirs(directory_name)
        reaction_rules_extraction(reaction_database_file, transformations, directory_name)
