from __future__ import annotations

__author__ = 'github.com/wardsimon'
__version__ = '0.0.1'

from easyCore.Objects.ObjectClasses import BaseObj, Parameter
from typing import ClassVar, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from easyCore.Utils.typing import iF


class Scattering(BaseObj):

    _name = 'scattering'
    _defaults = {
    }

    def __init__(self,
                 interface: Optional[iF] = None,
                 **kwargs):

        super().__init__(self._name, **kwargs)
        # Usually we don't want to fit these parameters
        # self.coherent.enabled = False

        self.interface = interface