======================
Spacegroups
======================

Constraints are a fundamental component in non-trivial fitting operations. They can also be used to affirm the minimum/maximum of a parameter or tie parameters together in a model.

Creating a spacegroup
-----------------------

Spacegroups can be created using a few class methods. The most common is by calling :class:`easyCrystallography.Components.SpaceGroup.SpaceGroup` directly, which will try to ascertain if the; Hermann-Mauguin name, Hall symbol, or international number was provided.

At the moment only known spacegroups can be used. The ability to describe custom spacegroups is planned for a future release.


Hermann-Mauguin notation
^^^^^^^^^^^^^^^^^^^^^^^^^
One of the extended Hermann-Mauguin symbols given in Table 4.3.2.1 of International Tables for Crystallography Vol. A (2002) can be used to generate a spacegroup.

.. code-block:: python

     from easyCrystallography.Components.SpaceGroup import SpaceGroup
     # Create a spacegroup from the Hermann-Mauguin symbol
     HM_symbol = 'C m c a'
     spacegroup = SpaceGroup(HM_symbol)
     print(spacegroup)
     # <Spacegroup: system: 'orthorhombic', number: 64, H-M: 'C m c a'>

IT Number
^^^^^^^^^^^^^^^^^^^^^^^^^
The number as assigned in International Tables for Crystallography Vol. A, specifying the proper affine class  of space groups
(crystallographic space-group type) to which the space group belongs.  This number defines the space-group type but not the coordinate system in which it is expressed. The method :meth:`easyCrystallography.Components.SpaceGroup.SpaceGroup.from_int_number` is used for generation.

.. code-block:: python

     from easyCrystallography.Components.SpaceGroup import SpaceGroup
     # Create a spacegroup from the IT number
     IT_number = 64
     spacegroup = SpaceGroup.from_int_number(IT_number)
     print(spacegroup)
     # <Spacegroup: system: 'orthorhombic', number: 64, H-M: 'C m c a'>

Hexagonal spacegroups can also be created by specifying the `hexagonal` flag.

.. code-block:: python

     from easyCrystallography.Components.SpaceGroup import SpaceGroup
     # Create a spacegroup from the IT number
     IT_number = 160
     spacegroup = SpaceGroup.from_int_number(IT_number, hexagonal=False)
     print(spacegroup)
     # <Spacegroup: system: 'trigonal', number: 160, H-M: 'R 3 m' setting: 'R'>
     spacegroup = SpaceGroup.from_int_number(IT_number, hexagonal=True)
     print(spacegroup)
     # <Spacegroup: system: 'trigonal', number: 160, H-M: 'R 3 m' setting: 'H'>

Hall Symbol
^^^^^^^^^^^^^^^^^^^^^^
Space-group symbol as described by Hall. This symbol gives the space-group setting explicitly. The spacegroup is found from these operators. So, arbitrary Hall symbols cannot be used to create a spacegroup.

.. code-block:: python

     from easyCrystallography.Components.SpaceGroup import SpaceGroup
     # Create a spacegroup from the Hall symbol
     Hall_symbol = '-C 2ac 2'
     spacegroup = SpaceGroup(HM_symbol)
     print(spacegroup)
     # <Spacegroup: system: 'orthorhombic', number: 64, H-M: 'C m c a'>

XYZ operators
^^^^^^^^^^^^^^^^^^^^^^
Often spacegroups are defined by a set `XYZ` operators. These operators are used to define the symmetry operations of the spacegroup. These can be used to create a spacegroup, but cannot be used to define arbitrarily groups. The method :meth:`easyCrystallography.Components.SpaceGroup.SpaceGroup.from_xyz_string` is used for generation.

.. code-block:: python

    from easyCrystallography.Components.SpaceGroup import SpaceGroup

    # Define the operators for the 'C m c a' spacegroup
    xyz = ['x,y,z',
            '-x,-y+1/2,z+1/2',
            'x,-y,-z',
            '-x,y+1/2,-z+1/2',
            '-x,-y,-z',
            'x,y+1/2,-z+1/2',
            '-x,y,z',
            'x,-y+1/2,z+1/2',
            'x+1/2,y+1/2,z',
            '-x+1/2,-y,z+1/2',
            'x+1/2,-y+1/2,-z',
            '-x+1/2,y,-z+1/2',
            '-x+1/2,-y+1/2,-z',
            'x+1/2,y,-z+1/2',
            '-x+1/2,y+1/2,z',
            'x+1/2,-y,z+1/2']
    xyz_string= ';'.join(xyz)
    spacegroup = SpaceGroup.from_xyz_string(xyz_string)
    print(spacegroup)
    # <Spacegroup: system: 'orthorhombic', number: 64, H-M: 'C m c a'>
    spacegroup = SpaceGroup.from_xyz_string(xyz)
    print(spacegroup)
    # <Spacegroup: system: 'orthorhombic', number: 64, H-M: 'C m c a'>

Symmetry matrices and operators
^^^^^^^^^^^^^^^^^^^^^^

Often spacegroups are defined by a set of rotation matrices and translation vectors. These can be used to create a spacegroup, but cannot be used to define arbitrarily groups. The method :meth:`easyCrystallography.Components.SpaceGroup.SpaceGroup.from_symMatrices`, :meth:`easyCrystallography.Components.SpaceGroup.SpaceGroup.from_generators` can be used for generation.

.. code-block:: python

    from easyCrystallography.Components.SpaceGroup import SpaceGroup

    # Define the operators for the 'P 2' spacegroup
    rot = [[[24, 0, 0], [0, 24, 0], [0, 0, 24]], [[-24, 0, 0], [0, 24, 0], [0, 0, -24]]]
    trans = [[0, 0, 0], [0, 0, 0]]

    spacegroup = SpaceGroup.from_symMatrices(rot, trans)
    print(spacegroup)
    # <Spacegroup: system: 'monoclinic', number: 3, H-M: 'P 1 2 1'>

    # You can also use generators
    rot = [[[24, 0, 0], [0, 24, 0], [0, 0, 24]], [[-24, 0, 0], [0, 24, 0], [0, 0, -24]]]
    trans = [[0, 0, 0]]
    spacegroup = SpaceGroup.from_generators(rot, trans)
    print(spacegroup)
    # <Spacegroup: system: 'monoclinic', number: 3, H-M: 'P 1 2 1'>

`easyCrystallography` uses it's own implementation of `Operations` (:class:`easyCrystallography.Symmetry.SymOp.SymmOp`) which can be used to describe arbitrary rotations and translations of points and objects. A spacegroup can also be formed from these operations. The following example uses rotations and translations to describe the spacegroup, however any of the class methods in the :class:`easyCrystallography.Symmetry.SymOp.SymmOp` can be used. The method :meth:`easyCrystallography.Components.SpaceGroup.SpaceGroup.from_symOps` is used for generation.

.. code-block:: python

    from easyCrystallography.Components.SpaceGroup import SpaceGroup
    from easyCrystallography.Symmetry.SymOp import SymmOp

    # Define the operators for the 'P 2' spacegroup
    rots = [[[24, 0, 0], [0, 24, 0], [0, 0, 24]], [[-24, 0, 0], [0, 24, 0], [0, 0, -24]]]
    trans = [[0, 0, 0], [0, 0, 0]]
    ops = []
    for rot, tran in zip(rots, trans):
        ops.append(SymmOp.from_rotation_and_translation(rot, tran))
    spacegroup = SpaceGroup.from_symOps(ops)
    print(spacegroup)
    # <Spacegroup: system: 'monoclinic', number: 3, H-M: 'P 1 2 1'>

The same can be used for `gemmi` operators and operator groups by using the :meth:`easyCrystallography.Components.SpaceGroup.SpaceGroup.from_gemmi_operations` class method.

Features
-----------------

Obtaining crystal symmetry information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The `SpaceGroup` class provides methods for obtaining crystal symmetry information. This includes the crystal system,
reference settings, unique identifiers and symmetry operations.

.. code-block:: python

    from easyCrystallography.Components.SpaceGroup import SpaceGroup

    # Create an interesting spacegroup.
    spacegroup_str = 'C m c a'
    spacegroup = SpaceGroup(spacegroup_str)

    # Obtain some information about the spacegroup.
    print(spacegroup.crystal_system)
    # 'orthorhombic'
    print(spacegroup.is_reference_setting)
    # True
    print(spacegroup.int_number)
    # 64
    print(spacegroup.hall_symbol)
    # '-C 2ac 2'


    # Obtain the symmetry operations.
    print(spacegroup.sym_xyz)
    # 'x, y, z;-x+1/2, -y, z+1/2;x, -y, -z;-x+1/2, y, -z+1/2;-x, -y, -z;x+1/2, y, -z+1/2;-x, y, z;x+1/2, -y, z+1/2;x+1/2, y+1/2, z;-x, -y+1/2, z+1/2;x+1/2, -y+1/2, -z;-x, y+1/2, -z+1/2;-x+1/2, -y+1/2, -z;x, y+1/2, -z+1/2;-x+1/2, y+1/2, z;x, -y+1/2, z+1/2'
    print(spacegroup.sym_ops)

Operating on a point
^^^^^^^^^^^^^^^^^^^^^^

We can use the operators in the spacegroup to transform points. This is useful for generating all atomic positions in a
structure. The following example uses the :meth:`easyCrystallography.Components.SpaceGroup.SpaceGroup.get_orbit` method.

.. code-block:: python

    from easyCore import np
    from easyCrystallography.Components.SpaceGroup import SpaceGroup

    # Create a spacegroup which has a few operations
    spacegroup_str = 'C m c a'
    spacegroup = SpaceGroup(spacegroup_str)

    # Create a point
    pt = np.array([0, 0.1, 0])
    pts = spacegroup.get_orbit(pt)
    print(pts)
    #[[ 0.   0.1  0. ]
    # [ 0.5 -0.1  0.5]
    # [ 0.  -0.1  0. ]
    # [ 0.5  0.1  0.5]
    # [ 0.5  0.6  0. ]
    # [ 0.   0.4  0.5]
    # [ 0.5  0.4  0. ]
    # [ 0.   0.6  0.5]]

Multiplicity for a point
^^^^^^^^^^^^^^^^^^^^^^^^
We can also use the spacegroup to determine the multiplicity of a point. This is useful for determining structure
factors and other properties. It is also required in some calculations. The following example uses the :meth:`easyCrystallography.Components.SpaceGroup.SpaceGroup.get_site_multiplicity` method.

.. code-block:: python

    from easyCore import np
    from easyCrystallography.Components.SpaceGroup import SpaceGroup

    # Create a spacegroup which has a few operations
    spacegroup_str = 'C m c a'
    spacegroup = SpaceGroup(spacegroup_str)

    # Create a point
    pt = np.array([0, 0.1, 0])
    mult = spacegroup.get_site_multiplicity(pt)
    print(mult)
    # 3


Symmetry Reference
--------------------

.. minigallery:: easyCore.Fitting.Constraints.NumericConstraint
    :add-heading: Examples using `SpaceGroup`