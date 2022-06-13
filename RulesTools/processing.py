from CGRtools import RDFRead, RDFWrite, ReactionContainer
from tqdm import tqdm
from collections import defaultdict
from pickle import dump


def apply_transformations(transformations: list, reaction: ReactionContainer) -> ReactionContainer:
    """
    Applies transformations to given reaction

    :param transformations: list of transformations
    :param reaction: the reaction to be transformed
    :return: result reaction
    """
    for transform in transformations:
        reaction = transform(reaction)
    return reaction


def apply_filters(reaction: ReactionContainer, reaction_filters):
    """

    :param reaction:
    :param reaction_filters:
    :return:
    """
    is_filtered = False
    for reaction_filter in reaction_filters:
        if reaction_filter(reaction):
            reaction.meta[reaction_filter.__class__.__name__] = 'True'
            is_filtered = True
    return reaction, is_filtered


def reaction_database_processing(reaction_database_file_name: str, transformations: list, filters: list, result_directory_name: str = './',
                                 error_reactions_file_name: str = 'error_reactions', rules_file_name: str = 'rules'):
    """
    Extracts reaction rules from the reaction database and applies other given transformations. Returns rules files in
    RDF and pickle formats and error reactions file in RDF format

    :param reaction_database_file_name: path to the reaction database
    :param transformations: list of transformations include rules extraction

    :param result_directory_name: result directory name
    :param error_reactions_file_name: error reactions file name
    :param rules_file_name: rules file name
    """
    # TODO опциональный выбор типа файла
    # TODO добавить распараллеливание?

    # def run_parallel(reactions, num_processes=38):
    #     with Pool(num_processes) as p:
    #         new_reactions = p.map(clean_reaction, reactions)
    #     return new_reactions


    unique_rules = defaultdict(lambda: [0, []])

    with RDFRead(reaction_database_file_name, indexable=True) as reactions:
        reactions.reset_index()

        for reaction in tqdm(reactions, total=len(reactions)):
            try:
                if filters:
                    reaction, is_filtered = apply_filters(reaction, filters)
                    print('filters good')
                    if is_filtered:
                        with RDFWrite(f'{result_directory_name}/{error_reactions_file_name}.rdf', append=True) as errors_file:
                            errors_file.write(reaction)
                if transformations:
                    rules = apply_transformations(transformations, reaction)
                    print('trans good')
            except Exception:
                with RDFWrite(f'{result_directory_name}/{error_reactions_file_name}.rdf', append=True) as errors_file:
                    reaction.meta['Error'] = 'True'
                    errors_file.write(reaction)
            else:
                for rule in rules:
                    unique_rules[rule][0] += 1
                    # unique_rules[rule][1].append(reaction.meta['Reaction_ID'])

    with RDFWrite(f'{result_directory_name}/{rules_file_name}.rdf') as rules_file:
        with open(f'{result_directory_name}/{rules_file_name}.pickle', 'ab') as rules_file2:
            for rule, info in unique_rules.items():
                rule.meta['number_of_reactions'] = info[0]
                # rule.meta['reactions_id'] = ' '.join(info[1])
                rule.clean2d()
                rules_file.write(rule)
                dump(rule, rules_file2)
