# -*- coding: utf-8 -*-


from pint import UnitRegistry, set_application_registry
ureg = UnitRegistry()
ureg.setup_matplotlib(True)
set_application_registry(ureg)