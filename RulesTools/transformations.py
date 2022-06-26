from CGRtools import MoleculeContainer, ReactionContainer, QueryContainer, Reactor
from typing import Literal, Union, List, Iterable
from itertools import islice, chain


class DeleteSmallMolecules:
    """Allows to remove "small" molecules from the reaction and save information about them"""

    def __init__(self, number_of_atoms: int = 6, small_molecules_to_meta: bool = False):
        """
        :param number_of_atoms: molecules with the number of atoms equal to number_of_atoms and below will be removed
        from the reaction
        :param small_molecules_to_meta: if True, deleted molecules are saved to meta
        """
        self.number_of_atoms = number_of_atoms
        self.small_molecules_to_meta = small_molecules_to_meta

    def __call__(self, reaction: ReactionContainer) -> ReactionContainer:
        """
        Removes molecules with the number of atoms equal to number_of_atoms or lower from the reaction

        :param reaction: a reaction object
        :return: the reaction without "small" molecules
        """
        new_reactants, small_reactants = self._split_molecules(reaction.reactants)
        new_products, small_products = self._split_molecules(reaction.products)

        new_reaction = ReactionContainer(new_reactants, new_products, reaction.reagents, reaction.meta)
        new_reaction.name = reaction.name

        if self.small_molecules_to_meta:
            united_small_reactants = self._molecules_to_molcontainer(small_reactants)
            new_reaction.meta['small_reactants'] = str(united_small_reactants)

            united_small_products = self._molecules_to_molcontainer(small_products)
            new_reaction.meta['small_products'] = str(united_small_products)

        return new_reaction

    def _split_molecules(self, molecules: Iterable) -> tuple:
        """
        Splits molecules according to the number of heavy atoms

        :param molecules:
        :return: "big" molecules and "small" molecules
        """
        big_molecules, small_molecules = [], []
        for molecule in molecules:
            if len(molecule) > self.number_of_atoms:
                big_molecules.append(molecule)
            else:
                small_molecules.append(molecule)

        return big_molecules, small_molecules

    @staticmethod
    def _molecules_to_molcontainer(molecules: Iterable) -> MoleculeContainer:
        """
        Unites molecules into one molecular container

        :param molecules: set of molecules
        :return: one molecular graph
        """
        molcontainer = MoleculeContainer()
        for molecule in molecules:
            molcontainer = molcontainer.union(molecule)

        return molcontainer


class RebalanceReaction:
    """Allows to rebalance the reaction by assembling CGR and then decomposing it"""

    def __call__(self, reaction: ReactionContainer) -> ReactionContainer:
        """
        Rebalances the reaction by assembling CGR and then decomposing it. Works for all reactions for which the correct
        CGR can be assembled

        :param reaction: a reaction object
        :return: a rebalanced reaction
        """
        cgr = ~reaction
        reactants, products = ~cgr
        rebalanced_reaction = ReactionContainer(reactants.split(), products.split(), reaction.reagents, reaction.meta)
        rebalanced_reaction.name = reaction.name
        return rebalanced_reaction


class ReverseReaction:
    """Allows to reverse reaction"""

    def __call__(self, reaction: ReactionContainer) -> ReactionContainer:
        """
        Reverses given reaction

        :param reaction: a reaction object
        :return: the reversed reaction
        """
        reversed_reaction = ReactionContainer(reaction.products, reaction.reactants, reaction.reagents, reaction.meta)
        reversed_reaction.name = reaction.name

        return reversed_reaction


class RemoveReagents:
    """Allows to remove reagents from reactants and products"""

    def __init__(self, keep_reagents: bool = True, reagents_max_size: int = 7):
        """
        :param keep_reagents: if True, the reagents are written to ReactionContainer
        :param reagents_max_size: max size of molecules that are called reagents, bigger are deleted
        """
        self.keep_reagents = keep_reagents
        self.regents_max_size = reagents_max_size

    def __call__(self, reaction: ReactionContainer) -> ReactionContainer:
        """
        Removes reagents (not changed molecules or molecules not involved in the reaction) from reactants and products

        :param reaction: a reaction object
        :return: cleaned reaction
        """
        new_reactants = list(reaction.reactants)
        new_products = list(reaction.products)
        new_reagents = list(reaction.reagents)

        not_changed_molecules = set(new_reactants).intersection(new_products)

        if not_changed_molecules:
            new_reactants = [molecule for molecule in new_reactants if molecule not in not_changed_molecules]
            new_products = [molecule for molecule in new_products if molecule not in not_changed_molecules]
            new_reagents = new_reagents.extend(not_changed_molecules)

        cgr = ~reaction
        center_atoms = set(cgr.center_atoms)

        for n, molecule in enumerate(new_reactants.copy()):
            if center_atoms.isdisjoint(molecule):
                new_reactants.pop(n)
                new_reagents.append(molecule)

        for n, molecule in enumerate(new_products.copy()):
            if center_atoms.isdisjoint(molecule):
                new_products.pop(n)
                new_reagents.append(molecule)

        if self.keep_reagents:
            new_reagents = [molecule for molecule in new_reagents if len(molecule) <= self.regents_max_size]
        else:
            new_reagents = []

        new_reaction = ReactionContainer(new_reactants, new_products, new_reagents, reaction.meta)
        new_reaction.name = reaction.name

        return new_reaction


class AddStandardize:
    """Allows to add reaction standardization"""

    def __call__(self, reaction: ReactionContainer) -> ReactionContainer:
        """
        Apply CGRtools standardize method

        :param reaction: a reaction object
        :return: standardized reaction
        """
        reaction.standardize()
        return reaction


class AddCanonicalize:
    """Allows to add reaction canonicalization"""

    def __call__(self, reaction: ReactionContainer) -> ReactionContainer:
        """
        Apply CGRtools canonicalize method

        :param reaction: a reaction object
        :return: canonicalized reaction
        """
        reaction.canonicalize()
        return reaction


class CreateRule:
    """Allows to create reaction rule of various types"""

    def __init__(self, rules_from_multistage_reaction: bool = False, environment_atoms_number: int = 1,
                 rule_with_functional_groups: bool = False,
                 functional_groups_list: List[Union[MoleculeContainer, QueryContainer]] = None,
                 include_rings: bool = False, add_leaving_and_coming: bool = False,
                 keep_reagents: bool = True, keep_meta: bool = True, as_query: bool = True,
                 keep_atom_info: Literal['none', 'reaction_center', 'all'] = 'reaction_center',
                 clean_info: Union[frozenset[str], str] = frozenset(
                     {'neighbors', 'hybridization', 'implicit_hydrogens', 'ring_sizes'}),
                 check_in_reactor: bool = False):
        """
        :param rules_from_multistage_reaction: if True, then it extracts all reaction rules of a given type separately
        from a multistep reaction
        :param environment_atoms_number: the reaction rule with the reaction center extended to the number of atoms
        equal environment_atoms_number is extracted
        :param rule_with_functional_groups: if True, then the reaction rule is extracted, the extended reaction center
        includes the functional groups adjacent to the reaction center
        :param functional_groups_list: list of functional groups contained as MoleculeContainer or QueryContainer
        :param include_rings: if True, it additionally extracts all rings adjacent to the reaction center
        :param add_leaving_and_coming: if True, it additionally adds leaving and coming groups to the rule
        :param as_query: if True, then returns a reaction rule with query containers and preserves or deletes specified
        structural properties of atoms such as: neighbors, hybridization, implicit hydrogens, ring sizes
        :param keep_atom_info: if 'all', then information about all atoms is saved, if 'reaction_center', then
        information is saved only about reaction center atoms, if 'none', then no information is saved
        :param clean_info: a frozenset of atom info types
        :param check_in_reactor: if True, the rule is applied to the reaction reactants and if the reaction products are
        not obtained, returns ValueError
        """
        assert set(clean_info) & {'neighbors', 'hybridization', 'implicit_hydrogens', 'ring_sizes'}

        if rule_with_functional_groups and not isinstance(functional_groups_list, Iterable):
            raise TypeError('Invalid type of functional group list')

        self.rules_from_multistage_reaction = rules_from_multistage_reaction
        self.environment_atoms_number = environment_atoms_number
        self.rule_with_functional_groups = rule_with_functional_groups
        self.functional_groups_list = functional_groups_list
        self.include_rings = include_rings
        self.add_leaving_and_coming = add_leaving_and_coming
        self.keep_reagents = keep_reagents
        self.keep_meta = keep_meta
        self.as_query = as_query
        self.keep_atom_info = keep_atom_info
        self.clean_info = clean_info
        self.check_in_reactor = check_in_reactor

    def __call__(self, reaction: ReactionContainer) -> List[ReactionContainer]:
        """
        Creates reaction rule from the reaction

        :param reaction: a reaction object
        :return: A reaction rule
        """
        if self.rules_from_multistage_reaction:
            reaction_rules = set()
            for single_reaction in islice(reaction.enumerate_centers(), 15):
                reaction_rule = self._rule_extraction(single_reaction)
                reaction_rules.add(reaction_rule)
            return list(reaction_rules)
        return [self._rule_extraction(reaction)]

    def _rule_extraction(self, reaction: ReactionContainer) -> ReactionContainer:
        """
        Creates reaction rule from the reaction

        :param reaction: a reaction object
        :return: A reaction rule
        """
        cgr = ~reaction
        center_atoms_numbers = set(cgr.center_atoms)
        rule_atoms_numbers = center_atoms_numbers.copy()

        if self.environment_atoms_number:
            reduced_cgr = cgr.augmented_substructure(center_atoms_numbers, deep=self.environment_atoms_number)
            rule_atoms_numbers = rule_atoms_numbers | set(reduced_cgr)

        if self.rule_with_functional_groups:
            for molecule in reaction.molecules():
                for functional_group in self.functional_groups_list:

                    # mapping of the functional group across the molecule
                    for mapping in functional_group.get_mapping(molecule):
                        functional_group.remap(mapping)

                        # check for intersection of the functional group with the reaction centre
                        if set(functional_group.atoms_numbers) & center_atoms_numbers:
                            rule_atoms_numbers |= set(functional_group.atoms_numbers)

                        remapping = {value: key for key, value in mapping.items()}
                        functional_group.remap(remapping)

        if self.include_rings:
            for ring in cgr.connected_rings:
                if set(ring) & center_atoms_numbers:
                    rule_atoms_numbers |= set(ring)

        if self.add_leaving_and_coming:
            reactants_atoms = {atom for reactant in reaction.reactants for atom in reactant}
            products_atoms = {atom for product in reaction.products for atom in product}
            leaving_atoms = reactants_atoms - products_atoms
            coming_atoms = products_atoms - reactants_atoms
            rule_atoms_numbers |= leaving_atoms | coming_atoms

        rule_reactants = [reactant.substructure(rule_atoms_numbers.intersection(reactant.atoms_numbers)) for reactant in
                          reaction.reactants if rule_atoms_numbers.intersection(reactant.atoms_numbers)]
        rule_products = [product.substructure(rule_atoms_numbers.intersection(product.atoms_numbers)) for product in
                         reaction.products if rule_atoms_numbers.intersection(product.atoms_numbers)]

        if self.as_query:
            rule_reactants = self._clean_rule_molecules(rule_reactants, reaction.reactants, center_atoms_numbers)
            rule_products = self._clean_rule_molecules(rule_products, reaction.products, center_atoms_numbers)

        if self.keep_reagents and self.as_query:
            rule_reagents = [reagent.substructure(reagent, as_query=True) for reagent in reaction.reagents]
        elif self.keep_reagents:
            rule_reagents = reaction.reagents
        else:
            rule_reagents = []

        if self.keep_meta:
            rule_meta = reaction.meta
        else:
            rule_meta = []

        reaction_rule = ReactionContainer(rule_reactants, rule_products, rule_reagents, rule_meta)
        reaction_rule.name = reaction.name
        reaction_rule.flush_cache()

        if self.check_in_reactor and self.as_query:
            reactor = Reactor(reaction_rule)
            reaction_products = list(reaction.products)
            reaction_products.sort()
            for result_reaction in reactor(reaction.reactants):
                result_products = []
                for result_product in result_reaction.products:
                    result_products.extend(result_product.split())
                result_products.sort()
                if result_products == reaction_products:
                    return reaction_rule
            raise ValueError('The reaction rule incorrectly generates products from reactants in the Reactor')

        return reaction_rule

    def _clean_rule_molecules(self, rule_molecules: (tuple, list), reaction_molecules: (tuple, list),
                              center_atoms: set) -> list:
        """
        For each rule molecule in the rule molecules list, it creates a query molecule and then removes specified
        information about atoms

        :param rule_molecules: a list of rule molecules
        :param reaction_molecules: a list of reaction molecules
        :param center_atoms: the atoms that are in the reaction center
        """
        cleaned_rule_molecules = []

        for rule_molecule in rule_molecules:
            for reaction_molecule in reaction_molecules:
                if set(rule_molecule.atoms_numbers).issubset(reaction_molecule.atoms_numbers):

                    query_reaction_molecule = reaction_molecule.substructure(reaction_molecule, as_query=True)
                    query_rule_molecule = query_reaction_molecule.substructure(rule_molecule)

                    if self.keep_atom_info == 'reaction_center':
                        for atom_number in set(rule_molecule.atoms_numbers) - set(center_atoms):
                            query_rule_molecule = self._clean_query_atom(query_rule_molecule, atom_number)
                    elif self.keep_atom_info == 'none':
                        for atom_number in rule_molecule.atoms_numbers:
                            query_rule_molecule = self._clean_query_atom(query_rule_molecule, atom_number)

                    cleaned_rule_molecules.append(query_rule_molecule)
                    break

        return cleaned_rule_molecules

    def _clean_query_atom(self, query_molecule: QueryContainer, atom_number: int) -> QueryContainer:
        """
        Removes the specified information about atom with atom_number from the query molecule

        :param query_molecule: the query molecule
        :param atom_number: the number of the atom to be modified
        :return: The query molecule with the atom information removed
        """
        for info_type in self.clean_info:
            if info_type == 'neighbors':
                query_molecule.atom(atom_number).neighbors = None
            if info_type == 'hybridization':
                query_molecule.atom(atom_number).hybridization = None
            if info_type == 'implicit_hydrogens':
                query_molecule.atom(atom_number).implicit_hydrogens = None
            if info_type == 'ring_sizes':
                query_molecule.atom(atom_number).ring_sizes = None

        return query_molecule


class CreateNotQuery:
    """Allows to create not query reaction (reaction rule) from query"""

    def __call__(self, rule: ReactionContainer) -> ReactionContainer:
        """
        Creates not query reaction (reaction rule) from query

        :param rule: reaction (reaction rule) with query molecules
        :return: reaction (reaction rule) with not query molecules
        """
        not_query_reactants = self.create_not_query_molecules(rule.reactants)
        not_query_products = self.create_not_query_molecules(rule.products)
        not_query_reagents = self.create_not_query_molecules(rule.reagents)
        not_query_rule = ReactionContainer(not_query_reactants, not_query_products, not_query_reagents, rule.meta)
        not_query_rule.name = rule.name
        return not_query_rule

    @staticmethod
    def create_not_query_molecules(query_molecules: Iterable) -> list[MoleculeContainer]:
        """
        Creates not query molecules from queries

        :param query_molecules: set of query molecules
        :return: set of not query molecules
        """
        not_query_molecules = []
        for molecule in query_molecules:
            not_query_molecule = MoleculeContainer()
            for atom_number, atom_type in molecule.atoms():
                not_query_molecule.add_atom(atom_type.atomic_number, atom_number, charge=atom_type.charge,
                                            is_radical=atom_type.is_radical)
            for atom_number1, atom_number2, bond_type in molecule.bonds():
                not_query_molecule.add_bond(atom_number1, atom_number2, int(bond_type))
            not_query_molecules.append(not_query_molecule)
        return not_query_molecules
