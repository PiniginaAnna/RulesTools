from CGRtools.containers import ReactionContainer, MoleculeContainer
from StructureFingerprint import MorganFingerprint
import numpy as np
from typing import Iterable
from CGRtools.containers.bonds import DynamicBond
from CGRtools.exceptions import InvalidAromaticRing


class CGRConnectedComponents:
    def __call__(self, reaction: ReactionContainer) -> bool:
        """
        Returns True if CGR contains unrelated components (without reagents), else False

        :param reaction: input reaction
        :return: True or False
        """
        tmp_reaction = ReactionContainer(reaction.reactants, reaction.products)  # нужно убрать реагенты
        cgr = ~tmp_reaction
        if cgr.connected_components_count > 1:
            return True
        else:
            return False


# TODO оно может убирать нормальные молекулы, надо ли
class IsStrangeCarbons:
    def __call__(self, reaction: ReactionContainer) -> bool:
        """
        Checks for the presence of methane or C molecules with only one type of bond (not aromatic) in the set of molecules

        :param reaction: input reaction
        :return: True or False
        """
        free_carbons = False
        for molecule in reaction.reactants + reaction.products:
            atoms_types = list(set(a.atomic_symbol for _, a in molecule.atoms()))  # atoms types in molecule
            if len(atoms_types) == 1:
                if atoms_types[0] == 'C':
                    if len(molecule) == 1:  # methane
                        free_carbons = True
                        break
                    else:
                        bond_types = list(set(int(b) for _, _, b in molecule.bonds()))
                        if len(bond_types) == 1:
                            if bond_types[0] != 4:
                                free_carbons = True  # C molecules with only one type of bond (not aromatic)
                                break
        return free_carbons


# TODO спросить в чем разница между ними
class IsCompeteReaction:
    def __call__(self, reaction: ReactionContainer) -> bool:
        """
        Returns True if there are compete reactions, else False

        :param reaction: input reaction
        :return: True or False
        """
        is_compete = False
        for mol in reaction.reactants + reaction.products:
            if len(mol) > 6:
                for reagent in reaction.reagents:
                    if len(reagent) > 6:
                        try:
                            clique_size = len(next(mol.get_mcs_mapping(reagent, limit=100)))
                            tanimoto = clique_size / (len(mol) + len(reagent) - clique_size)
                            if tanimoto > 0.6:
                                is_compete = True
                                break
                        except StopIteration:
                            continue
        return is_compete


def tanimoto_kernel(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Calculate Tanimoto between each elements of array x and y.
    Parameters
    ----------
    x : 2D array
        Array of features.
    y : 2D array
        Array of features.
    Note
    ----
    Features in arrays x and y should be equal and in same order.
    Returns
    -------
    array : 2D array
        Pairwise Tanimoto coefficients.
    """
    x_dot = np.dot(x, y.T)

    x2 = (x**2).sum(axis=1)
    y2 = (y**2).sum(axis=1)

    len_x2 = len(x2)
    len_y2 = len(y2)

    result = x_dot / (np.array([x2] * len_y2).T + np.array([y2] * len_x2) - x_dot)
    result[np.isnan(result)] = 0

    return result


class IsCompeteReactionMorgan:
    def __call__(self, reaction: ReactionContainer) -> bool:
        """
        Returns True if there are compete reactions, else False

        :param reaction: input reaction
        :return: True or False
        """
        mf = MorganFingerprint()
        is_compete = False
        for mol in reaction.reactants + reaction.products:
            if len(mol) > 6:
                for reagent in reaction.reagents:
                    if len(reagent) > 6:
                        molf = mf.transform([mol])
                        reagentf = mf.transform([reagent])
                        similarity = tanimoto_kernel(molf, reagentf)[0][0]
                        if len(mol) > 14 and len(reagent) > 14:
                            if similarity > 0.4:
                                is_compete = True
                                break
                        else:
                            if similarity > 0.1:
                                is_compete = True
                                break
        return is_compete


class CheckDynamicBondsNumber:
    def __init__(self, min_bonds_number: int = 1, max_bonds_number: int = 6):
        self.min_bonds_number = min_bonds_number
        self.max_bonds_number = max_bonds_number

    def __call__(self, reaction: ReactionContainer, ):
        cgr = ~reaction
        if self.min_bonds_number <= len(cgr.center_bonds) <= self.max_bonds_number:
            return False
        return True


class CheckSmallMolecules:
    def __init__(self, limit: int = 6):
        self.limit = limit

    def __call__(self, reaction: ReactionContainer):
        # проверка на маленькость
        if len(reaction.reactants) == 1 and self.is_only_small_molecules(reaction.reactants):  # только один маленький реактант
            return True
        elif len(reaction.products) == 1 and self.is_only_small_molecules(reaction.products):  # только один маленький продукт
            return True
        elif self.is_only_small_molecules(reaction.reactants) and self.is_only_small_molecules(reaction.products):  # только маленькие молекулы и в реактантах, и в продуктах
            return True
        return False

    def is_only_small_molecules(self, molecules: Iterable[MoleculeContainer]) -> bool:
        """
        Returns True if there are only small molecules in input, else False

        :param molecules: set of molecules
        :return: True or False
        """
        only_small_mols = True
        for molecule in molecules:
            if len(molecule) > self.limit:
                only_small_mols = False
                break
        return only_small_mols


class IsNoReaction:
    def __call__(self, reaction: ReactionContainer):
        # проверка на отсутствие реакции
        cgr = ~reaction
        if not cgr.center_atoms and not cgr.center_bonds:
            return True
        return False


class IsMultiCenter:
    def __call__(self, reaction: ReactionContainer):
        # проверка на мультицентровость
        cgr = ~reaction
        if len(cgr.centers_list) > 1:
            return True
        return False


class CheckRings:
    def __call__(self, reaction: ReactionContainer):
        # проверка валентности
        reaction.kekule()
        if reaction.check_valence():
            return True
        else:
            reaction.thiele()
            r_rings, r_arom_rings = self._calc_rings(reaction.reactants)
            p_rings, p_arom_rings = self._calc_rings(reaction.products)
            # не совпадает число ароматических колец в реактантах и продуктах
            if r_arom_rings != p_arom_rings:
                return True
            # не совпадает число колец в реактантах и продуктах
            elif r_rings != p_rings:
                return True
            else:
                return False

    @staticmethod
    def _calc_rings(molecules: Iterable) -> tuple[int, int]:
        """
        Calculates number of all rings and number of aromatic rings in molecules

        :param molecules: set of molecules
        :return: number of all rings and number of aromatic rings in molecules
        """
        rings, arom_rings = 0, 0
        for mol in molecules:
            rings += mol.rings_count
            arom_rings += len(mol.aromatic_rings)
        return rings, arom_rings


# TODO какая гибридизация у С? полный кгр? образование только одинарной с-с?
# TODO странный тест на карбены - проверка на наличие гетероатомов возле разорвавшейся с-н и образовавшейся с-с
# TODO проверка на конденсацию из-за с-н кислотности неправильная
class IsWrongCHBreaking:
    def __call__(self, reaction: ReactionContainer) -> bool:
        """
        Returns True if there are C-C bonds formation from breaking C-H exclude other than C-H at the alpha-position of the carbonyl group, to exclude
        condensation reactions and reactions with carbenes, else False

        :param reaction: input reaction
        :return: True or False
        """
        try:
            reaction.kekule()
            if not reaction.check_valence():
                reaction.thiele()
                reaction.explicify_hydrogens()
                cgr = ~reaction
                reduced_cgr = cgr.augmented_substructure(cgr.center_atoms, deep=1)
                if self.is_wrong_c_h_breaking(reduced_cgr):
                    reaction.meta['c_h_breaking'] = 'True'
            else:
                return True
        except InvalidAromaticRing:
            return True
        return False

    def is_wrong_c_h_breaking(self, cgr):
        """
        Returns True if there are C-C bonds formation from breaking C-H exclude other than C-H at the alpha-position of the carbonyl group, to exclude
        condensation reactions and reactions with carbenes, else False

        :param cgr: CGR with explicified hydrogens
        :return: True or False
        """
        neigh_carbon_id = None
        for cgr_atom_id, neighbours in cgr._bonds.items():
            for neighbour_id, bond in neighbours.items():
                if bond == DynamicBond(1, None):
                    center_atom = cgr.atom(cgr_atom_id)
                    neighbour = cgr.atom(neighbour_id)
                    if center_atom.atomic_symbol == 'C' and neighbour.atomic_symbol == 'H':  # разрывающаяся одинарная c-h связь

                        for oth_neighbour_id, oth_bond in neighbours.items():
                            if oth_neighbour_id != neighbour_id:
                                oth_neighbour = cgr.atom(oth_neighbour_id)
                                if oth_neighbour.atomic_symbol == 'C' and oth_bond == DynamicBond(None, 1):  # образовалась одинарная с-с
                                    # return True

                                    neigh_carbon_id = oth_neighbour_id

                                if oth_bond.order and oth_bond.p_order:  # not None orders
                                    condensation = oth_neighbour.atomic_symbol not in ('C', 'H')
                                    condensation = condensation and oth_bond.order == oth_bond.p_order
                                    condensation = condensation and oth_bond.order > 1
                                    if condensation:
                                        return False

                    if neigh_carbon_id:
                        ultimate_neighbours = list(cgr._bonds[neigh_carbon_id].keys()) + list(neighbours.keys())
                        for neigh_id in ultimate_neighbours:  # test on carbens?
                            if cgr.atom(neigh_id).atomic_symbol not in ('C', 'H'):
                                return False
                        return True
        return False


class IsCCsp3Breaking:
    def __call__(self, reaction: ReactionContainer) -> bool:
        """
        Returns True if there are C(sp3)-C(sp3) bonds breaking, else False

        :param reaction: input reaction
        :return: True or False
        """
        cgr = ~reaction
        reaction_center = cgr.augmented_substructure(cgr.center_atoms, deep=0)
        for atom_id, neighbour_id, bond in reaction_center.bonds():
            atom = reaction_center.atom(atom_id)
            neighbour = reaction_center.atom(neighbour_id)

            is_bond_broken = bond.order is not None and bond.p_order is None
            is_atoms_carbons = atom.atomic_symbol == 'C' and neighbour.atomic_symbol == 'C'
            is_atoms_sp3 = atom.hybridization == 1 and neighbour.hybridization == 1

            if is_bond_broken and is_atoms_carbons and is_atoms_sp3:
                return True
        return False


# TODO порядок связи в разрывающемся цикле?
class IsCCRingBreaking:
    def __call__(self, reaction: ReactionContainer) -> bool:
        """
        Returns True if there are ring C-C bonds breaking (5, 6, 7 atoms in rings), else False

        :param reaction: input reaction
        :return: True or False
        """
        cgr = ~reaction

        reactants_center_atoms = {}
        reactants_rings = ()
        for mol in reaction.reactants:
            reactants_rings += mol.sssr
            for n, atom in mol.atoms():
                if n in cgr.center_atoms:
                    reactants_center_atoms[n] = atom

        reaction_center = cgr.augmented_substructure(cgr.center_atoms, deep=0)
        for atom_id, neighbour_id, bond in reaction_center.bonds():
            atom = reactants_center_atoms[atom_id]
            neighbour = reactants_center_atoms[neighbour_id]

            is_bond_broken = bond.order is not None and bond.p_order is None
            is_atoms_carbons = atom.atomic_symbol == 'C' and neighbour.atomic_symbol == 'C'
            is_atoms_in_ring = bool({5, 6, 7}.intersection(atom.ring_sizes)) and \
                               bool({5, 6, 7}.intersection(neighbour.ring_sizes)) and \
                               any(atom_id in ring and neighbour_id in ring for ring in reactants_rings)

            if is_bond_broken and is_atoms_carbons and is_atoms_in_ring:
                return True
        return False
















# TODO может канониколайз?
# def clean_reaction(reaction: ReactionContainer, small_mol_size: int = 6) -> tuple[str, ReactionContainer]:
#     """
#
#     :param reaction:
#     :param small_mol_size:
#     :return:
#     """
#     try:
#         tmp_reaction = reaction.copy()
#         tmp_reaction.standardize()  # стандартизация групп
#
#         # удаление реагентов
#         tmp_reaction = remove_unused_molecules(tmp_reaction)
#         tmp_reaction = remove_reagents(tmp_reaction)
#
#         # проверка на маленькость
#         if len(tmp_reaction.reactants) == 1 and is_only_small_molecules(tmp_reaction.reactants, small_mol_size):  # один маленький реактант
#             return 'small_mol', reaction
#         elif len(tmp_reaction.products) == 1 and is_only_small_molecules(tmp_reaction.products, small_mol_size):  # один маленький продукт
#             return 'small_mol', reaction
#         elif is_only_small_molecules(tmp_reaction.reactants, small_mol_size) \
#                 and is_only_small_molecules(tmp_reaction.products, small_mol_size):  # только маленькие молекулы и в реактантах, и в продуктах
#             return 'small_mol', reaction
#
#         # проверка на конкурирующие продукты
#         # if is_compete_reaction(tmp_reaction.reactants, tmp_reaction.products, tmp_reaction.reagents):
#         if is_compete_reaction_morgan(tmp_reaction.reactants, tmp_reaction.products, tmp_reaction.reagents):
#             return 'compete', reaction
#
#         # проверка на несвязанные компоненты в КГР
#         if cgr_connected_components(tmp_reaction):
#             return 'cgr_connected', reaction
#
#         tmp2_reaction = ReactionContainer(tmp_reaction.reactants, tmp_reaction.products)
#
#         # проверка на отсутствие реакции
#         cgr = ~tmp2_reaction
#         if not cgr.center_atoms and not cgr.center_bonds:
#             return 'no_reaction', reaction
#
#         # проверка на мультицентровость
#         if len(cgr.centers_list) > 1:
#             return 'multi_center', reaction
#
#         # проверка на размер реакционного центра
#         reaction_center_size = len(cgr.center_atoms) + len(cgr.center_bonds)
#         reaction.meta['reaction_center_size'] = reaction_center_size
#         tmp_reaction.meta['reaction_center_size'] = reaction_center_size
#
#         if reaction_center_size <= 2 or reaction_center_size >= 15:
#             return 'reaction_center_size', reaction
#
#         # ребалансировка реакции
#         rebalanced_reaction = rebalance_reaction(tmp_reaction, cgr)
#
#         # проверка на "странные" углероды
#         if is_strange_carbons(list(rebalanced_reaction.reactants) + list(rebalanced_reaction.products)):
#             return 'carbons', tmp_reaction
#
#         # проверка валентности
#         rebalanced_reaction.kekule()
#         if rebalanced_reaction.check_valence():
#             return 'valence_mistake', tmp_reaction
#         else:
#             rebalanced_reaction.thiele()
#             rebalanced_reaction.clean2d()
#             r_rings, r_arom_rings = _calc_rings(rebalanced_reaction.reactants)
#             p_rings, p_arom_rings = _calc_rings(rebalanced_reaction.products)
#             # не совпадает число ароматических колец в реактантах и продуктах
#             if r_arom_rings != p_arom_rings:
#                 return 'not_equal_arom_rings', rebalanced_reaction
#             # не совпадает число колец в реактантах и продуктах
#             elif r_rings != p_rings:
#                 return 'not_equal_rings', rebalanced_reaction
#             else:
#                 return 'clean', rebalanced_reaction  # супер чистая реакция
#
#     except InvalidAromaticRing:
#         return 'invalid_arom_ring', reaction


# TODO зачем складывать атомы со связями? и зачем записывать в мету
# # проверка на размер реакционного центра
# reaction_center_size = len(cgr.center_atoms) + len(cgr.center_bonds)
# reaction.meta['reaction_center_size'] = reaction_center_size
# tmp_reaction.meta['reaction_center_size'] = reaction_center_size
#
# if reaction_center_size <= 2 or reaction_center_size >= 15:
#     return 'reaction_center_size', reaction
