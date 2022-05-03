from CGRtools import RDFRead, RDFWrite, ReactionContainer
from tqdm import tqdm
from collections import defaultdict
from pickle import dump


def _apply_transformations(transformations: list, reaction: ReactionContainer) -> ReactionContainer:
    """
    Applies transformations to given reaction

    :param transformations: list of transformations
    :param reaction: the reaction to be transformed
    :return: result reaction
    """
    for transform in transformations:
        reaction = transform(reaction)
    return reaction


def reaction_rules_extraction(reaction_database_file_name: str, transformations: list, directory_name: str = './',
                              error_reactions_file_name: str = 'error_reactions', rules_file_name: str = 'rules'):
    """
    Extracts reaction rules from the reaction database and applies other given transformations. Returns rules files in
    RDF and pickle formats and error reactions file in RDF format

    :param reaction_database_file_name: path to the reaction database
    :param transformations: list of transformations include reaction extraction
    :param directory_name: result directory name
    :param error_reactions_file_name: error reactions file name
    :param rules_file_name: rules file name
    """
    unique_rules = defaultdict(lambda: [0, []])

    with RDFRead(reaction_database_file_name, indexable=True) as reactions:
        reactions.reset_index()

        for reaction in tqdm(reactions, total=len(reactions)):
            try:
                rules = _apply_transformations(transformations, reaction)
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
