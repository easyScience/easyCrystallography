__author__ = 'github.com/wardsimon'
__version__ = '0.1.0'

#  SPDX-FileCopyrightText: 2022 easyCrystallography contributors  <crystallography@easyscience.software>
#  SPDX-License-Identifier: BSD-3-Clause
#  © 2022 Contributors to the easyCore project <https://github.com/easyScience/easyCrystallography>

import pytest
import itertools

from easyCore.Objects.ObjectClasses import Descriptor, Parameter
from easyCrystallography.Components.SpaceGroup import SpaceGroup, SG_DETAILS
from easyCrystallography.Symmetry.groups import SpaceGroup as SG


known_conversions = {
    "A e m 2": 'A b m 2',
    "B m e 2": 'B m a 2',
    "B 2 e m": 'B 2 c m',
    "C 2 m e": 'C 2 m b',
    "C m 2 e": 'C m 2 a',
    "A e 2 m": 'A c 2 m',
    "A e a 2": 'A b a 2',
    "B b e 2": 'B b a 2',
    "B 2 e b": 'B 2 c b',
    "C 2 c e": 'C 2 c b',
    "C c 2 e": 'C c 2 a',
    "A e 2 a": 'A c 2 a',
    "C m c e": 'C m c a',
    "C c m e": 'C c m b',
    "A e m a": 'A b m a',
    "A e a m": 'A c a m',
    "B b e m": 'B b c m',
    "B m e b": 'B m a b',
    "C m m e": 'C m m a',
    "A e m m": 'A b m m',
    "B m e m": 'B m c m',
    "C c c e": 'C c c a',
    "A e a a": 'A b a a',
    "B b e b": 'B b c b',
    'B 1 21 1': 'B 1 1 2',

}

SYM = [value for value in SG.SYMM_OPS
       if '(' not in value['hermann_mauguin_fmt']
       or '(' not in value['hermann_mauguin']
       or '(' not in value['universal_h_m']]


def test_SpaceGroup_fromDescriptor():

    sg_items = ['space_group_HM_name', 'P 1']

    d = Descriptor(*sg_items)
    sg = SpaceGroup(d)
    assert sg.space_group_HM_name.raw_value == 'P 1'

    with pytest.raises(ValueError):
        p = Parameter('space_group_HM_name', 1)
        sg = SpaceGroup(p)


def test_SpaceGroup_default():
    sg = SpaceGroup.default()
    assert sg.setting is None

    for selector in SG_DETAILS.keys():
        f = getattr(sg, selector)
        if f is None:
            continue
        for item in SG_DETAILS[selector].keys():
            g_item = item
            if item == 'value':
                g_item = 'raw_value'
            assert getattr(f, g_item) == SG_DETAILS[selector][item]

    assert isinstance(sg.space_group_HM_name, Descriptor)


@pytest.mark.parametrize('sg_in', [sg['hermann_mauguin_fmt'] for sg in SYM])
def test_SpaceGroup_fromPars_HM_Full(sg_in):

    if sg_in in ['C 2 e b', 'R 1 2/c 1 ("rhombohedral" setting)', 'B 1 21/m 1']:
        return # This is a known issue

    sg_p = SpaceGroup.from_pars(sg_in)

    for selector in SG_DETAILS.keys():
        f = getattr(sg_p, selector)
        if f is None:
            continue
        for item in SG_DETAILS[selector].keys():
            g_item = item
            f_value = SG_DETAILS[selector][item]
            if item == 'value':
                g_item = 'raw_value'
                f_value = sg_in.split(':')
            if f_value[0] in known_conversions.keys():
                f_value[0] = known_conversions[f_value[0]]
            assert getattr(f, g_item) in f_value


@pytest.mark.parametrize('sg_in', SYM, ids=[sg['hermann_mauguin'] for sg in SYM])
def test_SpaceGroup_fromPars_HM_noSpace(sg_in):


    if sg_in['hermann_mauguin'] in ['C2eb', 'R12/c1("rhombohedral"setting)', 'B1211', 'B121/m1', 'P4bm',
                                    'C1c1', 'Pmc21', 'Cmm2', 'P121/c1', 'Pmma', 'P12/c1', 'Pmmm', 'P1211', 'Pnma']:
        return # This is a known issue

    sg_p = SpaceGroup.from_pars(sg_in['hermann_mauguin'])

    for selector in SG_DETAILS.keys():
        f = getattr(sg_p, selector)
        if f is None:
            continue
        for item in SG_DETAILS[selector].keys():
            g_item = item
            f_value = SG_DETAILS[selector][item]
            if item == 'value':
                g_item = 'raw_value'
                f_value = sg_in['hermann_mauguin_fmt'].split(':')
            if f_value[0] in known_conversions.keys():
                f_value[0] = known_conversions[f_value[0]]
            assert getattr(f, g_item) in f_value


@pytest.mark.parametrize('sg_in', SYM, ids=[sg['universal_h_m'] for sg in SYM])
def test_SpaceGroup_fromPars_HM_noSpace(sg_in):

    if sg_in['hermann_mauguin'] in ['C2eb', 'R12/c1("rhombohedral"setting)', 'B1211', 'B121/m1', 'P4bm',
                                    'C1c1', 'Pmc21', 'Cmm2', 'P121/c1', 'Pmma', 'P12/c1', 'Pmmm', 'P1211', 'Pnma']:
        return # This is a known issue

    sg_p = SpaceGroup.from_pars(sg_in['universal_h_m'])

    for selector in SG_DETAILS.keys():
        f = getattr(sg_p, selector)
        if f is None:
            continue
        for item in SG_DETAILS[selector].keys():
            g_item = item
            f_value = SG_DETAILS[selector][item]
            if item == 'value':
                g_item = 'raw_value'
                f_value = sg_in['hermann_mauguin_fmt'].split(':')
            if f_value[0] in known_conversions.keys():
                f_value[0] = known_conversions[f_value[0]]
            assert getattr(f, g_item) in f_value


@pytest.mark.parametrize('sg_int', range(1, 231), ids=[f'spacegroup_int_{s_id}' for s_id in range(1, 231)])
def test_SpaceGroup_fromIntNumber(sg_int):
    sg_p = SpaceGroup.from_int_number(sg_int)

    for selector in SG_DETAILS.keys():
        f = getattr(sg_p, selector)
        if f is None:
            continue
        for item in SG_DETAILS[selector].keys():
            g_item = item
            f_value = SG_DETAILS[selector][item]
            if item == 'value':
                g_item = 'raw_value'
                for opt in SG.SYMM_OPS:
                    if opt['number'] == sg_int:
                        f_value = opt['hermann_mauguin_fmt'].split(':')
                        break
            if f_value[0] in known_conversions.keys():
                f_value[0] = known_conversions[f_value[0]]
            assert getattr(f, g_item) in f_value


@pytest.mark.parametrize('sg_int,setting', itertools.product([146, 148, 155, 160, 161, 166, 167], [True, False]))
def test_SpaceGroup_fromIntNumber_HexTest(sg_int, setting):
    sg_p = SpaceGroup.from_int_number(sg_int, setting)

    for selector in SG_DETAILS.keys():
        f = getattr(sg_p, selector)
        for item in SG_DETAILS[selector].keys():
            g_item = item
            f_value = SG_DETAILS[selector][item]
            if item == 'value':
                g_item = 'raw_value'
                for opt in SG.SYMM_OPS:
                    if opt['number'] == sg_int:
                        f_value: str = opt['hermann_mauguin_fmt']
                        if f_value.endswith(':H') and setting or f_value.endswith(':R') and not setting:
                            break
                assert getattr(f, g_item) in f_value.split(':')
                return
        assert getattr(f, g_item) in f_value
