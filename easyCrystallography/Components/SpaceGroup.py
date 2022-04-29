#  SPDX-FileCopyrightText: 2022 easyCrystallography contributors  <crystallography@easyscience.software>
#  SPDX-License-Identifier: BSD-3-Clause
#  © 2022 Contributors to the easyCore project <https://github.com/easyScience/easyCrystallography>
#
from __future__ import annotations

__author__ = "github.com/wardsimon"
__version__ = "0.1.0"

import gemmi
import numbers
from typing import ClassVar, Type, List, Tuple, Optional, TYPE_CHECKING, Union, Dict, NoReturn

from easyCrystallography.Symmetry.SymOp import SymmOp

from easyCore import np
from easyCore.Objects.ObjectClasses import BaseObj, Descriptor
from easyCore.Utils.io.star import StarEntry, StarSection, FakeCore, FakeItem, ItemHolder


SG_DETAILS = {
    "space_group_HM_name": {
        "name":        "hermann_mauguin",
        "description": "Hermann-Mauguin symbols given in Table 4.3.2.1 of International Tables for Crystallography "
                       "Vol. A (2002) or a Hermann-Mauguin symbol for a conventional or unconventional setting.",
        "url":         "https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ispace_group_name_H-M_alt.html",
        "value":       "P 1",
    },
    'setting':             {
        'name':        'coordinate-code',
        'description': 'A qualifier taken from the enumeration list identifying which setting in International Tables '
                       'for Crystallography Volume A (2002) (IT) is used.',
        'url':         'https://www.iucr.org/__data/iucr/cifdic_html/2/cif_sym.dic/Ispace_group'
                       '.IT_coordinate_system_code.html',
        'value':       '\x00',
    }
}

if TYPE_CHECKING:
    from easyCore.Utils.typing import iF
    from easyCrystallography.Components.Site import S
    import numpy.typing as npt
    T = Union[S, npt.ArrayLike]


class SpaceGroup(BaseObj):
    _space_group_HM_name: ClassVar[Descriptor]
    _setting: ClassVar[Descriptor]
    _CIF_SECTION_NAME: ClassVar[str] = "_space_group"
    _CIF_CONVERSIONS: ClassVar[List[Tuple[str, str]]] = [
        ("hermann_mauguin", "_name_H-M_ref"),
        ("setting_str", "_IT_coordinate_system_code"),
    ]

    def __init__(self,
                 space_group_HM_name: Optional[Descriptor, str] = None,
                 setting: Optional[Descriptor, str] = None,
                 interface: Optional[iF] = None):
        """
        Generate a spacegroup object from it's Hermann-Mauguin symbol and setting. The setting can be a part of the
        Hermann-Mauguin symbol.

        :param space_group_HM_name: Hermann-Mauguin symbol
        :param setting: Optional setting for the space group
        :param interface: Interface to the calculator
        """

        super(SpaceGroup, self).__init__(
            "space_group",
            _space_group_HM_name=Descriptor(**SG_DETAILS["space_group_HM_name"]),
            _setting=Descriptor(**SG_DETAILS["setting"]),
        )

        if space_group_HM_name:
            self._space_group_HM_name = space_group_HM_name
        if setting is not None:
            self._setting = setting
        self._operations = None
        self.__on_change(self._space_group_HM_name.raw_value, self._setting.raw_value, set_internal=True)
        self.interface = interface
        self._cell = None

    @classmethod
    def from_pars(cls, space_group_HM_name: str, setting: str = None, interface: Optional[iF] = None):
        """
        Generate a spacegroup object from it's Hermann-Mauguin symbol and setting.

        :param space_group_HM_name: Hermann-Mauguin symbol
        :param setting: Optional setting for the space group
        :param interface: Interface to the calculator
        """
        return cls(space_group_HM_name, setting, interface)

    @classmethod
    def default(cls, interface: Optional[iF] = None):
        """
        Generate a P1 spacegroup object.

        :param interface: Interface to the calculator
        """
        return cls(interface=interface)

    @classmethod
    def from_int_number(cls, int_number: int, hexagonal=True, interface: Optional[iF] = None):
        """
        Generate a spacegroup object from it's spacegroup number (1-231).
        :param int_number: spacegroup number
        :param hexagonal: Should a hexagonal setting be used?
        :param interface: Interface to the calculator
        """
        sg = gemmi.find_spacegroup_by_number(int(int_number))
        setting = None
        if int_number in [146, 148, 155, 160, 161, 166, 167]:
            if hexagonal:
                setting = 'H'
            else:
                setting = 'R'
        return cls(sg.hm, setting, interface=interface)

    @classmethod
    def from_symOps(cls, sym_ops: List[SymmOp], interface: Optional[iF] = None):
        """
        Create a space group from a list of easyCrystallography symmetry operations.

        :param sym_ops: List of easyCrystallography symmetry operations
        :param interface: Interface to the calculator
        """
        ops = [gemmi.Op(op.as_xyz_string()) for op in sym_ops]
        sg_data = gemmi.find_spacegroup_by_ops(ops)
        return cls(sg_data.hm, sg_data.ext, interface=interface)

    @classmethod
    def from_symMatrices(cls, rotations: List[np.ndarray], translations: List[np.ndarray], interface: Optional[iF] = None):
        """
        Create a space group from a lists of rotations and translations. The number of rotations and translations must
        be the same.

        :param rotations: Array of rotation matrices [n*[3x3]]
        :param translations: Array of translations [n*[1x3]]
        :param interface: Interface to the calculator
        """
        ops = []
        for rot, tran in zip(rotations, translations):
            ops.append(SymmOp.from_rotation_and_translation(rot, tran))
        return cls.from_symOps(ops, interface=interface)

    @classmethod
    def from_generators(cls, rotations: List[np.ndarray], translations: List[np.ndarray], interface: Optional[iF] = None):
        """
        Create a space group from a lists of rotations and translations. Each translation is applied to each rotation.

        :param rotations: Array of rotation matrices [n*[3x3]]
        :param translations: Array of translations [m*[1x3]]
        :param interface: Interface to the calculator
        """
        rots = len(translations)*rotations
        trans = [np.array(d) for d in np.array(translations).repeat(len(rotations), axis=0).tolist()]
        return cls.from_symMatrices(rots, trans, interface=interface)

    @classmethod
    def from_gemmi_operations(cls, operations: Union[gemmi.GroupOps, List[gemmi.Op]], interface: Optional[iF] = None):
        """
        Generate a space group from a list of gemmi operations or a gemmi group of operations.

        :param operations: Operations which define a space group
        :param interface: Interface to the calculator
        """
        ops = []
        if isinstance(operations, gemmi.GroupOps):
            ops = operations
        else:
            for op in operations:
                if isinstance(op, gemmi.Op):
                    ops.append(op)
                else:
                    raise TypeError("Operations must be of type gemmi.Op")
            ops = gemmi.GroupOps(ops)
        sg_data = gemmi.find_spacegroup_by_ops(ops)
        return cls(sg_data.hm, sg_data.ext, interface=interface)

    @classmethod
    def from_xyz_string(cls, xyz_string: Union[str, List[str]], interface: Optional[iF] = None):
        """
        Create a space group from a string or list of strings of the form 'x,y,z'.

        :param xyz_string: String defining space group operators
        :param interface: Interface to the calculator
        """
        if isinstance(xyz_string, str):
            xyz_string = xyz_string.split(';')
        ops = []
        for xyz in xyz_string:
            ops.append(SymmOp.from_xyz_string(xyz))
        return cls.from_symOps(ops, interface)

    def __on_change(self, new_spacegroup: Union[int, str],
                    new_setting: Optional = None,
                    set_internal: bool = True) -> Tuple[str, str]:
        """
        Internal function to update the space group. This function is called when the space group is changed. It checks
        the form of the imputs, generates reference data, and updates the internal data (if requested).

        :param new_spacegroup: New space group number or name
        :param new_setting: New space group setting
        :param set_internal: Should internal objects be updated
        """
        setting = "\x00"
        if isinstance(new_spacegroup, str):
            if ':' in new_spacegroup:
                new_spacegroup, new_setting = new_spacegroup.split(':')
            sg_data = gemmi.find_spacegroup_by_name(new_spacegroup)
            if sg_data is None:
                try:
                    sg_data = gemmi.find_spacegroup_by_ops(gemmi.symops_from_hall(new_spacegroup))
                except RuntimeError:
                    sg_data = None
        else:
            sg_data = gemmi.find_spacegroup_by_number(int(new_spacegroup))

        if sg_data is None:
            raise ValueError(f"Spacegroup \'{new_spacegroup}\' not found in database.")

        reference = sg_data.ext
        if new_setting is None or new_setting == "" or new_setting == "\x00":
            if reference != '\x00':
                setting = reference
        else:
            if new_setting != reference:
                sg_data = gemmi.find_spacegroup_by_name(sg_data.hm + ':' + new_setting)
                setting = sg_data.ext
            else:
                setting = new_setting
        if set_internal:
            self._sg_data = sg_data
            self._space_group_HM_name.value = sg_data.hm
            self._setting.value = setting
            self._operations = [SymmOp.from_rotation_and_translation(op.rot, op.tran) for op in sg_data.operations()]
        return sg_data.hm, setting

    @property
    def setting(self) -> Optional[Descriptor]:
        """
        Space group setting. If the space group does not have a setting, this will be None.

        :return: Space group setting
        """
        setting_str = self._setting.raw_value
        if setting_str == "\x00":
            return None  # no setting
        return self._setting

    @setting.setter
    def setting(self, new_setting: str) -> NoReturn:
        """
        Set the space group setting.

        :param new_setting: Space group setting
        """
        _, setting = self.__on_change(self._space_group_HM_name.raw_value, new_setting, set_internal=True)

    @property
    def setting_str(self) -> str:
        """
        Space group setting as a string. If the space group does not have a setting, this will be an empty string.

        :return: Space group setting as a string
        """
        if self.setting is None:
            return ""
        return self._setting.raw_value

    @property
    def space_group_HM_name(self) -> Descriptor:
        """
        Space group name as defined by a Hermann-Mauguin symbol

        :return: Space group name as easyCore Descriptor
        """
        return self._space_group_HM_name

    @space_group_HM_name.setter
    def space_group_HM_name(self, value: str) -> NoReturn:
        """
        Set the space group name as defined by a Hermann-Mauguin symbol

        :param value: Space group name as a string
        """
        _, _ = self.__on_change(value, set_internal=True)

    @property
    def hermann_mauguin(self) -> str:
        """
        Space group name as defined by a Hermann-Mauguin symbol

        :return: Space group name as a string
        """
        return self._space_group_HM_name.raw_value

    @property
    def hall_symbol(self) -> str:
        """
        Hall symbol of the space group

        :return: Hall symbol of the space group
        """
        return self._sg_data.hall

    @property
    def int_number(self) -> int:
        """
        International number of the space group

        :return: International number of the space group
        """
        return self._sg_data.number

    @int_number.setter
    def int_number(self, new_it_number: int) -> NoReturn:
        """
        Set the spacegroup by its international number

        :param new_it_number: International number of the new space group
        """
        _, _ = self.__on_change(new_it_number, set_internal=True)

    @property
    def crystal_system(self) -> str:
        """
        Which crystal system the space group belongs to

        :return: Crystal system of the space group
        """
        return self._sg_data.crystal_system_str()

    @property
    def symmetry_opts(self) -> List[SymmOp]:
        """
        List of symmetry operations of the space group

        :return: List of symmetry operations of the space group
        """
        return self._operations

    @property
    def symmetry_xyz(self) -> str:
        """
        Symmetry operations of the space group as a string

        :return: String of symmetry operations of the space group
        """
        return ';'.join([op.as_xyz_string() for op in self._operations])

    @property
    def is_reference_setting(self) -> bool:
        """
        Is the space group setting the reference setting?

        :return: Is the space group setting the reference setting?
        """
        return self._sg_data.is_reference_setting()

    def symmetry_matrices(self) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        """
        Get the rotational and translational matrices of the space group

        :return: Rotation and translation matrices
        """
        ops = self.symmetry_opts
        return [op.rotation_matrix.copy() for op in ops], [op.translation_vector.copy() for op in ops]

    def get_orbit(self, point: T, tol: float = 1e-5) -> np.ndarray:
        """
        Returns the orbit for a point.

        :param point: Point to get the orbit for
        :param tol: Tolerance for the orbit
        :return: Orbits of the point
        """
        orbit = []
        if not hasattr(point, '__iter__'):
            point = point.fract_coords
        point = np.array(point)
        for o in self.symmetry_opts:
            pp = np.array(o.operate(point))
            if not in_array_list(orbit, pp, tol=tol):
                orbit.append(pp)
        return np.array(orbit)

    def get_site_multiplicity(self, site: T, tol=1e-5) -> int:
        """
        Get the multiplicity of a given site

        :param site: Site to get the multiplicity of
        :param tol: Tolerance for the orbit
        :return: Multiplicity of the site
        """
        if not hasattr(site, '__iter__'):
            site = site.fract_coords
        site = np.array(site)
        multiplicity = 1
        for o in self.symmetry_opts:
            new_site = o.operate(site)
            if np.isclose(new_site, site, atol=tol).all():
                multiplicity += 1
        return multiplicity

    def to_star(self) -> StarEntry:
        """
        Convert the space group to a star entry

        :return: Star entry of the space group
        """
        if self.setting_str:
            s = FakeCore()
            item = FakeItem(self.hermann_mauguin)
            item.name = "space_group_HM_name"
            s._kwargs["space_group_HM_name"] = item
            item = FakeItem(self.setting_str)
            item.name = 'setting'
            s._kwargs['setting'] = item
            return StarSection(s)
        return StarEntry(self.space_group_HM_name)

    @classmethod
    def from_star(cls, in_string: str):
        """

        :param in_string:
        :type in_string:
        :return:
        :rtype:
        """
        return StarEntry.from_string(cls, in_string)

    def as_dict(self, skip: Optional[List[str]] = None) -> Dict[str, Union[str, bool, numbers.Number]]:
        d = super(SpaceGroup, self).as_dict(skip=skip)
        for item in self._CIF_CONVERSIONS:
            test = '_' + item[0]
            if test in d.keys():
                d[item[0]] = d.pop(test)
        return d

    @classmethod
    def from_cif_block(cls, block: gemmi.cif.Block):
        kwargs = {}
        for item in cls._CIF_CONVERSIONS:
            kwargs[item[0]] = block.find_pair_item(cls._CIF_SECTION_NAME + item[1])
        return cls(**kwargs)

    def to_cif_str(self) -> str:
        block = gemmi.cif.Block('temp')
        start_str = block.as_string()
        self.add_to_cif_block(block)
        final_str = block.as_string()
        return final_str[len(start_str):]

    def add_to_cif_block(self, block: gemmi.cif.Block) -> NoReturn:
        for item in self._CIF_CONVERSIONS:
            value = str(getattr(self, item[0]))
            if value:
                block.set_pair(self._CIF_SECTION_NAME + item[1], value)

    def __repr__(self) -> str:
        out_str = "<Spacegroup: system: '{:s}', number: {}, H-M: '{:s}'".format(
            self.crystal_system, self.int_number, self.hermann_mauguin
        )
        if self.setting_str:
            out_str = "{:s} setting: '{:s}'".format(out_str, self.setting_str)
        return out_str + ">"


def in_array_list(array_list, a, tol=1e-5) -> bool:
    """
    Extremely efficient nd-array comparison using numpy's broadcasting. This
    function checks if a particular array a, is present in a list of arrays.
    It works for arrays of any size, e.g., even matrix searches.

    Args:
        array_list ([array]): A list of arrays to compare to.
        a (array): The test array for comparison.
        tol (float): The tolerance. Defaults to 1e-5. If 0, an exact match is
            done.

    Returns:
        (bool)
    """
    if len(array_list) == 0:
        return False
    axes = tuple(range(1, a.ndim + 1))
    if not tol:
        return np.any(np.all(np.equal(array_list, a[None, :]), axes))
    return np.any(np.sum(np.abs(array_list - a[None, :]), axes) < tol)
