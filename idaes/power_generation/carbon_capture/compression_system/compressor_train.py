##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
import pyomo.environ as pyo
from pyomo.environ import (ConcreteModel,
                           SolverFactory,
                           TransformationFactory
                           )
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.power_generation.properties.natural_gas_PR import get_prop
import idaes.generic_models.properties.swco2 as swco2
from idaes.generic_models.unit_models import (
    Feed, Valve, ValveFunctionType, Heater)
from idaes.generic_models.unit_models.flash import Flash
# from compressor_PR import (CompressionStage, VaneDiffuserType, ImpellerType)
from idaes.power_generation.carbon_capture.compression_system.compressor_PR \
      import (CompressionStage, VaneDiffuserType, ImpellerType)

import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

from idaes.core.solvers import use_idaes_solver_configuration_defaults
import idaes


def main(m):

    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt["options"]["nlp_scaling_method"] = "user-scaling"
    # due to a lot of component mole fractions being on their lower bound of 0
    # bound push result in much longer solve times, so set it low.
    idaes.cfg.ipopt["options"]["bound_push"] = 1e-12

    m.fs = FlowsheetBlock(default={'dynamic': False})

    # Properties
    NG_config = get_prop(components=["H2O", 'CO2'], phases=["Vap"])
    NG_config_liq = get_prop(components=["H2O", 'CO2'], phases=["Vap", "Liq"])

    m.fs.NG_props = GenericParameterBlock(default=NG_config)
    m.fs.NG_props_liq = GenericParameterBlock(default=NG_config_liq)
    m.fs.props_co2 = swco2.SWCO2ParameterBlock(default={
                  "phase_presentation": swco2.PhaseType.G
              })

    # add feed
    m.fs.feed_air1 = Feed(default={"property_package": m.fs.NG_props_liq})

    # add valve
    m.fs.inlet_valve = Valve(
        default={"valve_function_callback": ValveFunctionType.linear,
                 "property_package": m.fs.NG_props,
                 "has_phase_equilibrium": True})

    m.fs.heater1 = Heater(default={"property_package": m.fs.NG_props_liq,
                                   "has_pressure_change": True,
                                   "has_phase_equilibrium": True})
    m.fs.heater2 = Heater(default={"property_package": m.fs.NG_props_liq,
                                   "has_pressure_change": True,
                                   "has_phase_equilibrium": True})
    m.fs.heater3 = Heater(default={"property_package": m.fs.NG_props_liq,
                                   "has_pressure_change": True,
                                   "has_phase_equilibrium": True})
    m.fs.heater4 = Heater(default={"property_package": m.fs.NG_props_liq,
                                   "has_pressure_change": True,
                                   "has_phase_equilibrium": True})
    m.fs.heater5 = Heater(default={"property_package": m.fs.NG_props_liq,
                                   "has_pressure_change": True,
                                   "has_phase_equilibrium": True})
    m.fs.heater6 = Heater(default={"property_package": m.fs.NG_props_liq,
                                   "has_pressure_change": True,
                                   "has_phase_equilibrium": True})
    m.fs.heater7 = Heater(default={"property_package": m.fs.props_co2,
                                   "has_pressure_change": True,
                                   "has_phase_equilibrium": True})

    m.fs.compstage1 = CompressionStage(
      default={"property_package": m.fs.NG_props,
               "impeller_type": ImpellerType.open_impeller,
               "vane_diffuser_type": VaneDiffuserType.vane_diffuser})
    m.fs.compstage2 = CompressionStage(
      default={"property_package": m.fs.NG_props,
               "impeller_type": ImpellerType.open_impeller,
               "vane_diffuser_type": VaneDiffuserType.vane_diffuser})
    m.fs.compstage3 = CompressionStage(
      default={"property_package": m.fs.NG_props,
               "impeller_type": ImpellerType.open_impeller,
               "vane_diffuser_type": VaneDiffuserType.vane_diffuser})
    m.fs.compstage4 = CompressionStage(
      default={"property_package": m.fs.NG_props,
               "impeller_type": ImpellerType.open_impeller,
               "vane_diffuser_type": VaneDiffuserType.vane_diffuser})
    m.fs.compstage5 = CompressionStage(
        default={"property_package": m.fs.NG_props,
                 "impeller_type": ImpellerType.open_impeller,
                 "vane_diffuser_type": VaneDiffuserType.vane_diffuser})
    m.fs.compstage6 = CompressionStage(
        default={"property_package": m.fs.NG_props,
                 "impeller_type": ImpellerType.open_impeller,
                 "vane_diffuser_type": VaneDiffuserType.vane_diffuser})
    m.fs.compstage7 = CompressionStage(
        default={"property_package": m.fs.props_co2,
                 "impeller_type": ImpellerType.open_impeller,
                 "vane_diffuser_type": VaneDiffuserType.vane_diffuser})
    m.fs.compstage8 = CompressionStage(
        default={"property_package": m.fs.props_co2,
                 "impeller_type": ImpellerType.open_impeller,
                 "vane_diffuser_type": VaneDiffuserType.vane_diffuser})

    m.fs.flash1 = Flash(default={"property_package": m.fs.NG_props_liq,
                                 'has_pressure_change': True})
    m.fs.flash2 = Flash(default={"property_package": m.fs.NG_props_liq,
                                 'has_pressure_change': True})
    m.fs.flash3 = Flash(default={"property_package": m.fs.NG_props_liq,
                                 'has_pressure_change': True})
    m.fs.flash4 = Flash(default={"property_package": m.fs.NG_props_liq,
                                 'has_pressure_change': True})
    m.fs.flash5 = Flash(default={"property_package": m.fs.NG_props_liq,
                                 'has_pressure_change': True})
    m.fs.flash6 = Flash(default={"property_package": m.fs.NG_props_liq,
                                 'has_pressure_change': True})

    # Arc
    m.fs.feedFG = Arc(
      source=m.fs.feed_air1.outlet, destination=m.fs.inlet_valve.inlet)
    m.fs.s1 = Arc(source=m.fs.inlet_valve.outlet,
                  destination=m.fs.compstage1.inlet)
    m.fs.s2 = Arc(
      source=m.fs.compstage1.outlet, destination=m.fs.heater1.inlet)
    m.fs.s3 = Arc(
      source=m.fs.heater1.outlet, destination=m.fs.flash1.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling flow
    iscale.set_scaling_factor(
      m.fs.compstage1.control_volume.properties_in[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage1.control_volume.properties_out[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage2.control_volume.properties_in[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage2.control_volume.properties_out[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage3.control_volume.properties_in[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage3.control_volume.properties_out[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage4.control_volume.properties_in[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage4.control_volume.properties_out[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage5.control_volume.properties_in[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage5.control_volume.properties_out[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage6.control_volume.properties_in[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.compstage6.control_volume.properties_out[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.heater5.control_volume.properties_in[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.heater5.control_volume.properties_out[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.flash5.control_volume.properties_in[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.flash5.control_volume.properties_out[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.heater6.control_volume.properties_in[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.heater6.control_volume.properties_out[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.flash6.control_volume.properties_in[0].flow_mol, 1e-6)
    iscale.set_scaling_factor(
      m.fs.flash6.control_volume.properties_out[0].flow_mol, 1e-6)

    # call scaling
    iscale.calculate_scaling_factors(m)

    # deactivate unit to avoid convergence issue
    # run simulation and initialize each compressor stage with its intercooler
    m.fs.compstage2.deactivate()
    m.fs.heater2.deactivate()
    m.fs.flash2.deactivate()
    m.fs.compstage3.deactivate()
    m.fs.heater3.deactivate()
    m.fs.flash3.deactivate()
    m.fs.compstage4.deactivate()
    m.fs.heater4.deactivate()
    m.fs.flash4.deactivate()
    m.fs.compstage5.deactivate()
    m.fs.heater5.deactivate()
    m.fs.flash5.deactivate()
    m.fs.compstage6.deactivate()
    m.fs.heater6.deactivate()
    m.fs.flash6.deactivate()
    m.fs.compstage7.deactivate()
    m.fs.heater7.deactivate()
    m.fs.compstage8.deactivate()

    # flowsheet
    # Variables
    # feed inlet
    p1 = 1.15 * 1e5  # Pa
    t1 = 40.0113 + 273.15  # K
    f1 = 1579.798470673201  # mol/s
    # m.fs.feed_air1.flow_mol[:] = f1*1e5
    m.fs.feed_air1.flow_mol.fix(f1)
    m.fs.feed_air1.temperature.fix(t1)
    m.fs.feed_air1.pressure.fix(p1)
    m.fs.feed_air1.mole_frac_comp[:, "H2O"].fix(0.0601355)
    m.fs.feed_air1.mole_frac_comp[:, "CO2"].fix(0.939864)
    m.fs.feed_air1.initialize()

    propagate_state(m.fs.feedFG)
    deltaP = (1.15-1.050328097156878)*1e5
    valve_opening = 48.8387/100
    cv = f1/pyo.sqrt(deltaP)/valve_opening

    m.fs.inlet_valve.Cv.fix(cv)
    m.fs.inlet_valve.valve_opening.fix(valve_opening)
    m.fs.inlet_valve.initialize()

    propagate_state(m.fs.s1)
    # fix compressor specification
    m.fs.compstage1.speed_of_sound.fix(279.7815969751224)
    m.fs.compstage1.heat_capacity_ratio.fix(1.2772)
    m.fs.compstage1.U2.fix(315.3)
    m.fs.compstage1.r2.fix(0.67654)
    m.fs.compstage1.z_s.fix(0.9952)
    m.fs.compstage1.z_d1.fix(0.97373)
    m.fs.compstage1.efficiency_mech.fix(0.97)
    m.fs.compstage1.eff_drive.fix(1.0)

    m.fs.compstage1.initialize(outlvl=idaeslog.INFO_HIGH)

    m.fs.compstage1_coeff_a = pyo.Var(
        m.fs.time, doc="coefficient a")  # initialize=-71.881,
    m.fs.compstage1_coeff_b = pyo.Var(
        m.fs.time, doc="coefficient b")  # initialize=10.6219,
    m.fs.compstage1_coeff_c = pyo.Var(
        m.fs.time, doc="coefficient c")  # initialize=0.805509,
    m.fs.compstage1_Ang = pyo.Var(m.fs.time, doc="IGV")  # initialize=62.8746,

    @m.fs.Constraint(m.fs.time)
    def coeff_a_eqn(b, t):
        return b.compstage1_coeff_a[t] == -38.46 * pyo.exp(
            -0.005098 * b.compstage1_Ang[t]) + (- 0.1043) * pyo.exp(
                0.09821 * b.compstage1_Ang[t])

    @m.fs.Constraint(m.fs.time)
    def coeff_b_eqn(b, t):
        return b.compstage1_coeff_b[t] == 6.863 * pyo.exp(
              -0.01329 * b.compstage1_Ang[t]) \
            + 0.02018 * pyo.exp(0.09443 * b.compstage1_Ang[t])

    @m.fs.Constraint(m.fs.time)
    def coeff_c_eqn(b, t):
        return b.compstage1_coeff_c[t] == -0.00891 * pyo.exp(
            0.06757 * b.compstage1_Ang[t]) + 0.9661 * pyo.exp(
                0.006228 * b.compstage1_Ang[t])

    @m.fs.Constraint(
        m.fs.time, doc="Correlation between psi_s and psi_3")
    def psi_s_stage1_eqn(b, t):
        return b.compstage1.psi_s[t] == b.compstage1_coeff_a[t] \
            * b.compstage1.psi_3[t]**2 + b.compstage1_coeff_b[
                t] * b.compstage1.psi_3[t] + b.compstage1_coeff_c[t]

    m.fs.compstage1.ratioP.fix(2.375530641906031)

    propagate_state(m.fs.s2)
    deltaP_1 = -0.009249261595310397*1e5 - 0.1157197593400905*1e5  # Pa
    heat_duty_1 = -6891953.465092183 * 1.076  # J/s
    m.fs.heater1.deltaP.fix(deltaP_1)
    m.fs.heater1.heat_duty.fix(heat_duty_1)
    m.fs.heater1.initialize()

    propagate_state(m.fs.s3)
    m.fs.flash1.deltaP.fix(0)
    m.fs.flash1.heat_duty.fix(0)
    m.fs.flash1.initialize()

    # remove bounds of variables, avoid convergence issue
    strip_bounds = pyo.TransformationFactory("contrib.strip_var_bounds")
    strip_bounds.apply_to(m, reversible=False)

    print('degrees of freedoms =', degrees_of_freedom(m.fs))
    solver = SolverFactory('ipopt')
    solver.options = {"max_iter": 500}
    solver.solve(m.fs, tee=True)

    m.fs.compstage2.activate()
    m.fs.V1 = Arc(
      source=m.fs.flash1.vap_outlet, destination=m.fs.compstage2.inlet)
    m.fs.heater2.activate()
    m.fs.s4 = Arc(
      source=m.fs.compstage2.outlet, destination=m.fs.heater2.inlet)
    m.fs.flash2.activate()
    m.fs.s5 = Arc(
      source=m.fs.heater2.outlet, destination=m.fs.flash2.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    propagate_state(m.fs.V1)
    m.fs.compstage2.speed_of_sound.fix(274.9242056090937)
    m.fs.compstage2.heat_capacity_ratio.fix(1.2772)
    m.fs.compstage2.U2.fix(315.3)
    m.fs.compstage2.r2.fix(0.507405)
    m.fs.compstage2.z_s.fix(0.9952)
    m.fs.compstage2.z_d1.fix(0.97373)
    m.fs.compstage2.efficiency_mech.fix(0.97)
    m.fs.compstage2.eff_drive.fix(1.0)

    m.fs.compstage2.initialize(outlvl=idaeslog.INFO_HIGH)
    m.fs.compstage2.ratioP.fix(2.424235107713805)

    m.fs.compstage2_coeff_a = pyo.Var(
        initialize=-36.1, doc="coefficient a")
    m.fs.compstage2_coeff_b = pyo.Var(
        initialize=5.8604, doc="coefficient b")
    m.fs.compstage2_coeff_c = pyo.Var(
        initialize=0.9737, doc="coefficient c")
    m.fs.compstage2_coeff_a.fix()
    m.fs.compstage2_coeff_b.fix()

    @m.fs.Constraint(m.fs.time, doc="Correlation between psi_s and psi_3")
    def psi_s_stage2_eqn(b, t):
        return b.compstage2.psi_s[t] == \
            b.compstage2_coeff_a * b.compstage2.psi_3[
                t]**2 + b.compstage2_coeff_b * b.compstage2.psi_3[
                    t] + b.compstage2_coeff_c

    propagate_state(m.fs.s4)
    m.fs.heater2.heat_duty.fix(0)
    m.fs.heater2.deltaP.fix(0)

    # Initialize the model.
    m.fs.heater2.initialize()
    propagate_state(m.fs.s5)
    deltaP_2 = -0.008968346342640359 * 1e5 - 0.08872504908530339 * 1e5  # Pa
    duty_2 = -5171342.002714329 * 1.09  # J/s
    m.fs.flash2.deltaP.fix(deltaP_2)
    m.fs.flash2.heat_duty.fix(duty_2)
    m.fs.flash2.initialize()

    solver.solve(m.fs, tee=True)

    m.fs.compstage3.activate()
    m.fs.V2 = Arc(
      source=m.fs.flash2.vap_outlet, destination=m.fs.compstage3.inlet)
    m.fs.heater3.activate()
    m.fs.s6 = Arc(
      source=m.fs.compstage3.outlet, destination=m.fs.heater3.inlet)
    m.fs.flash3.activate()
    m.fs.s7 = Arc(
      source=m.fs.heater3.outlet, destination=m.fs.flash3.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    propagate_state(m.fs.V2)
    m.fs.compstage3.speed_of_sound.fix(271.9618977806346)
    m.fs.compstage3.heat_capacity_ratio.fix(1.2766)
    m.fs.compstage3.U2.fix(280.4)
    m.fs.compstage3.r2.fix(0.3038799999999998)
    m.fs.compstage3.z_s.fix(0.97373)
    m.fs.compstage3.z_d1.fix(0.88949)
    m.fs.compstage3.efficiency_mech.fix(0.97)
    m.fs.compstage3.eff_drive.fix(1.0)

    m.fs.compstage3.initialize(outlvl=idaeslog.INFO_HIGH)

    m.fs.compstage3_coeff_a = pyo.Var(
        initialize=-33.47, doc="coefficient a")
    m.fs.compstage3_coeff_b = pyo.Var(
        initialize=6.619, doc="coefficient b")
    m.fs.compstage3_coeff_c = pyo.Var(
        initialize=0.9557, doc="coefficient c")
    m.fs.compstage3_coeff_a.fix()
    m.fs.compstage3_coeff_b.fix()
    # m.fs.compstage3_coeff_c.fix()

    @m.fs.Constraint(m.fs.time, doc="Correlation between psi_s and psi_3")
    def psi_s_stage3_eqn(b, t):
        return b.compstage3.psi_s[t] == \
            b.compstage3_coeff_a * b.compstage3.psi_3[
                t]**2 + b.compstage3_coeff_b * b.compstage3.psi_3[
                    t] + b.compstage3_coeff_c

    m.fs.compstage3.ratioP.fix(2.079358277332422)

    propagate_state(m.fs.s6)
    m.fs.heater3.heat_duty.fix(0)
    m.fs.heater3.deltaP.fix(0)

    # Initialize the model.
    m.fs.heater3.initialize()

    propagate_state(m.fs.s7)
    deltaP_3 = -0.008924593036713899 * 1e5 - 0.1132042852130902 * 1e5  # Pa
    duty_3 = -4034541.444294178 * 1.052  # J/s
    m.fs.flash3.deltaP.fix(deltaP_3)
    m.fs.flash3.heat_duty.fix(duty_3)
    m.fs.flash3.initialize()

    solver.solve(m.fs, tee=True)

    m.fs.compstage4.activate()
    m.fs.V3 = Arc(
      source=m.fs.flash3.vap_outlet, destination=m.fs.compstage4.inlet)
    m.fs.heater4.activate()
    m.fs.flash4.activate()
    m.fs.s8 = Arc(
      source=m.fs.compstage4.outlet, destination=m.fs.heater4.inlet)
    m.fs.s9 = Arc(
      source=m.fs.heater4.outlet, destination=m.fs.flash4.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    propagate_state(m.fs.V3)
    m.fs.compstage4.speed_of_sound.fix(267.5297492192697)
    m.fs.compstage4.heat_capacity_ratio.fix(1.2766)
    m.fs.compstage4.U2.fix(280.4)
    m.fs.compstage4.r2.fix(0.2431039999999999)
    m.fs.compstage4.z_s.fix(0.97373)
    m.fs.compstage4.z_d1.fix(0.88949)
    m.fs.compstage4.efficiency_mech.fix(0.97)
    m.fs.compstage4.eff_drive.fix(1.0)

    m.fs.compstage4.initialize(outlvl=idaeslog.INFO_HIGH)

    m.fs.compstage4_coeff_a = pyo.Var(
        initialize=-30.985, doc="coefficient a")
    m.fs.compstage4_coeff_b = pyo.Var(
        initialize=4.5087, doc="coefficient b")
    m.fs.compstage4_coeff_c = pyo.Var(
        initialize=0.9737, doc="coefficient c")
    m.fs.compstage4_coeff_a.fix()
    m.fs.compstage4_coeff_b.fix()
    # m.fs.compstage4_coeff_c.fix()

    @m.fs.Constraint(m.fs.time, doc="Correlation between psi_s and psi_3")
    def psi_s_stage4_eqn(b, t):
        return b.compstage4.psi_s[t] == \
            b.compstage4_coeff_a * b.compstage4.psi_3[
                t]**2 + b.compstage4_coeff_b * b.compstage4.psi_3[
                    t] + b.compstage4_coeff_c

    m.fs.compstage4.ratioP.fix(1.942059278698412)

    propagate_state(m.fs.s8)
    m.fs.heater4.heat_duty.fix(0)
    m.fs.heater4.deltaP.fix(0)

    # Initialize the model.
    m.fs.heater4.initialize()

    propagate_state(m.fs.s9)
    deltaP_4 = -0.008841956712291546*1e5 - 0.1035281931558739*1e5   # Pa
    duty_4 = -3806787.210729675 * 1.05                              # J/s
    m.fs.flash4.deltaP.fix(deltaP_4)
    m.fs.flash4.heat_duty.fix(duty_4)
    m.fs.flash4.initialize()

    solver.solve(m.fs, tee=True)

    m.fs.compstage5.activate()

    m.fs.V4 = Arc(
      source=m.fs.flash4.vap_outlet, destination=m.fs.compstage5.inlet)
    m.fs.heater5.activate()
    m.fs.s10 = Arc(
      source=m.fs.compstage5.outlet, destination=m.fs.heater5.inlet)
    m.fs.flash5.activate()
    m.fs.s11 = Arc(
      source=m.fs.heater5.outlet, destination=m.fs.flash5.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.compstage5.speed_of_sound.fix(259.618461791984)
    m.fs.compstage5.heat_capacity_ratio.fix(1.275)
    m.fs.compstage5.U2.fix(228.6244982815907)
    m.fs.compstage5.r2.fix(0.1529078325394176)
    m.fs.compstage5.z_s.fix(0.88949)
    m.fs.compstage5.z_d1.fix(0.58755)
    m.fs.compstage5.efficiency_mech.fix(0.97)
    m.fs.compstage5.eff_drive.fix(1.0)

    m.fs.compstage5.initialize(outlvl=idaeslog.INFO_HIGH)

    m.fs.compstage5_coeff_a = pyo.Var(
        initialize=-30.93, doc="coefficient a")
    m.fs.compstage5_coeff_b = pyo.Var(
        initialize=6.857, doc="coefficient b")
    m.fs.compstage5_coeff_c = pyo.Var(
        initialize=0.9216999999999999, doc="coefficient c")
    m.fs.compstage5_coeff_a.fix()
    m.fs.compstage5_coeff_b.fix()

    @m.fs.Constraint(m.fs.time, doc="Correlation between psi_s and psi_3")
    def psi_s_stage5_eqn(b, t):
        return b.compstage5.psi_s[t] == \
            b.compstage5_coeff_a * b.compstage5.psi_3[t]**2 +\
            b.compstage5_coeff_b * b.compstage5.psi_3[t] + b.compstage5_coeff_c

    m.fs.compstage5.ratioP.fix(1.727832026596874)

    propagate_state(m.fs.s10)
    m.fs.heater5.heat_duty.fix(0)
    m.fs.heater5.deltaP.fix(0)

    # Initialize the model.
    m.fs.heater5.initialize()

    propagate_state(m.fs.s11)
    deltaP_5 = -0.009452313515632591*1e5 - 2.38525160398048*1e5  # Pa
    duty_5 = -3417625.676517817 * 1.05  # J/s

    m.fs.flash5.deltaP.fix(deltaP_5)
    m.fs.flash5.heat_duty.fix(duty_5)
    m.fs.flash5.initialize()

    solver.solve(m.fs, tee=True)

    m.fs.compstage6.activate()

    m.fs.V5 = Arc(
      source=m.fs.flash5.vap_outlet, destination=m.fs.compstage6.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.compstage6.speed_of_sound.fix(247.1753695900838)
    m.fs.compstage6.heat_capacity_ratio.fix(1.275)
    m.fs.compstage6.U2.fix(228.6244982815907)
    m.fs.compstage6.r2.fix(0.1223262660315341)
    m.fs.compstage6.z_s.fix(0.88949)
    m.fs.compstage6.z_d1.fix(0.58755)
    m.fs.compstage6.efficiency_mech.fix(0.97)
    m.fs.compstage6.eff_drive.fix(1.0)

    m.fs.compstage6.initialize(outlvl=idaeslog.INFO_HIGH)

    m.fs.compstage6_coeff_a = pyo.Var(
        initialize=-22.27, doc="coefficient a")
    m.fs.compstage6_coeff_b = pyo.Var(
        initialize=5.2785, doc="coefficient b")
    m.fs.compstage6_coeff_c = pyo.Var(
        initialize=0.7965999999999999, doc="coefficient c")
    m.fs.compstage6_coeff_a.fix()
    m.fs.compstage6_coeff_b.fix()

    @m.fs.Constraint(m.fs.time, doc="Correlation between psi_s and psi_3")
    def psi_s_stage6_eqn(b, t):
        return b.compstage6.psi_s[t] == \
            b.compstage6_coeff_a * b.compstage6.psi_3[t]**2 +\
            b.compstage6_coeff_b * b.compstage6.psi_3[t] + b.compstage6_coeff_c

    m.fs.compstage6.ratioP.fix(1.646794262021896)

    solver.solve(m.fs, tee=True)

    m.fs.heater6.activate()
    m.fs.s12 = Arc(
      source=m.fs.compstage6.outlet, destination=m.fs.heater6.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    propagate_state(m.fs.s12)

    m.fs.heater6.heat_duty.fix(0)
    m.fs.heater6.deltaP.fix(0)

    # Initialize the model.
    m.fs.heater6.initialize()

    solver.solve(m.fs, tee=True)

    m.fs.flash6.activate()
    m.fs.s13 = Arc(
      source=m.fs.heater6.outlet, destination=m.fs.flash6.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    propagate_state(m.fs.s13)
    deltaP_6 = -0.009615469329299659*1e5 - 0.146127384189479*1e5  # Pa
    duty_6 = -4620568.766244856 * 1.05  # J/s
    m.fs.flash6.deltaP.fix(deltaP_6)
    m.fs.flash6.heat_duty.fix(duty_6)
    m.fs.flash6.initialize()

    solver.solve(m.fs, tee=True)

    # outlet of flash 6 will be connected to the surrogate TEG model.
    # because that IDAES model is not available yet,
    # assuming that the same feed inlet from Aspen for compressor stage 7
    m.fs.compstage7.activate()

    p7 = 57.72116303818244*1e5  # Pa
    t7 = 305.7688277980252  # K
    f7 = 1452.145808833169  # mol/s
    h7 = swco2.htpx(T=t7*pyo.units.K, P=p7*pyo.units.Pa)  # J/mol
    m.fs.compstage7.inlet.flow_mol.fix(f7)
    m.fs.compstage7.inlet.enth_mol.fix(h7)
    m.fs.compstage7.inlet.pressure.fix(p7)
    m.fs.compstage7.speed_of_sound.fix(222.4103311128723)
    m.fs.compstage7.heat_capacity_ratio.fix(1.275)  # 1.275
    m.fs.compstage7.U2.fix(192.4317606614612)
    m.fs.compstage7.r2.fix(0.07499999999999999)
    m.fs.compstage7.z_s.fix(0.58755)
    m.fs.compstage7.z_d1.fix(0.58755)
    m.fs.compstage7.efficiency_mech.fix(0.97)
    m.fs.compstage7.eff_drive.fix(1.0)

    m.fs.compstage7.initialize(outlvl=idaeslog.INFO_HIGH)

    m.fs.compstage7_coeff_a = pyo.Var(
        initialize=-16.04, doc="coefficient a")
    m.fs.compstage7_coeff_b = pyo.Var(
        initialize=6.038, doc="coefficient b")
    m.fs.compstage7_coeff_c = pyo.Var(
        initialize=0.6715899999999999, doc="coefficient c")
    m.fs.compstage7_coeff_a.fix()
    m.fs.compstage7_coeff_b.fix()

    @m.fs.Constraint(m.fs.time, doc="Correlation between psi_s and psi_3")
    def psi_s_stage7_eqn(b, t):
        return b.compstage7.psi_s[t] == \
            b.compstage7_coeff_a * b.compstage7.psi_3[t]**2 +\
            b.compstage7_coeff_b * b.compstage7.psi_3[t] + b.compstage7_coeff_c

    m.fs.compstage7.ratioP.fix(1.76868)

    m.fs.heater7.activate()
    m.fs.compstage8.activate()
    m.fs.s15 = Arc(
      source=m.fs.compstage7.outlet, destination=m.fs.heater7.inlet)
    m.fs.s16 = Arc(
      source=m.fs.heater7.outlet, destination=m.fs.compstage8.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    propagate_state(m.fs.s15)
    deltaP_7 = 0  # Pa
    heat_duty_7 = -1358006.853155364  # J/s
    m.fs.heater7.deltaP.fix(deltaP_7)
    m.fs.heater7.heat_duty.fix(heat_duty_7)
    m.fs.heater7.initialize()

    propagate_state(m.fs.s16)

    m.fs.compstage8.speed_of_sound.fix(246.0586429400477)
    m.fs.compstage8.heat_capacity_ratio.fix(1.275)
    m.fs.compstage8.U2.fix(192.4317606614612)
    m.fs.compstage8.r2.fix(0.06375)
    m.fs.compstage8.z_s.fix(0.58755)
    m.fs.compstage8.z_d1.fix(0.58755)
    m.fs.compstage8.efficiency_mech.fix(0.97)
    m.fs.compstage8.eff_drive.fix(1.0)

    m.fs.compstage8.initialize(outlvl=idaeslog.INFO_HIGH)

    m.fs.compstage8_coeff_a = pyo.Var(
        initialize=-16.04, doc="coefficient a")
    m.fs.compstage8_coeff_b = pyo.Var(
        initialize=5.3485, doc="coefficient b")
    m.fs.compstage8_coeff_c = pyo.Var(
        initialize=0.5515, doc="coefficient c")
    m.fs.compstage8_coeff_a.fix()
    m.fs.compstage8_coeff_b.fix()

    @m.fs.Constraint(m.fs.time, doc="Correlation between psi_s and psi_3")
    def psi_s_stage8_eqn(b, t):
        return b.compstage8.psi_s[t] == \
            b.compstage8_coeff_a * b.compstage8.psi_3[t]**2 +\
            b.compstage8_coeff_b * b.compstage8.psi_3[t] + b.compstage8_coeff_c

    m.fs.compstage8.ratioP.fix(1.4966)

    # strip_bounds = pyo.TransformationFactory("contrib.strip_var_bounds")
    # strip_bounds.apply_to(m, reversible=False)

    print('degrees of freedoms =', degrees_of_freedom(m.fs))
    if degrees_of_freedom(m.fs) == 0:
        solver = SolverFactory('ipopt')
        solver.options = {"max_iter": 500}
        solver.solve(m.fs, tee=True)
    else:
        raise Exception("Degrees of freedom is not 0.")

    print('Tin of stage 7 =', pyo.value(
        m.fs.compstage7.control_volume.properties_in[0].temperature))
    print('Tout of stage 7=', pyo.value(
        m.fs.compstage7.control_volume.properties_out[0].temperature))
    print('Pout of stage 7=', pyo.value(
        m.fs.compstage7.control_volume.properties_out[0].pressure))
    print('Tout of cooler 7 =', pyo.value(
        m.fs.heater7.control_volume.properties_out[0].temperature))
    print('Tout of product =', pyo.value(
        m.fs.compstage8.control_volume.properties_out[0].temperature))
    print('Pout of product =', pyo.value(
        m.fs.compstage8.control_volume.properties_out[0].pressure))

    return m


if __name__ == "__main__":
    m = ConcreteModel()
    main(m)
