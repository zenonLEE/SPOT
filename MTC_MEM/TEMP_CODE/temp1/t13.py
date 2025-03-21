import numpy as np
from thermo import ChemicalConstantsPackage

from insulator import Insulation
from prop_calculator import VLEThermo

a = VLEThermo(["CO2", "H2", "Methanol", "H2O", "carbon monoxide"])

gas = [0.00250396, 0.00782662, 0.00167475, 0.00108385, 0.00015737]
gas = np.array(gas)
# print(gas / np.sum(gas))
Tc = 323
Din = 0.02
Dd = 0.002
r = 0.95
Do = Din + Dd * 2
print(a.phi(503, 70, gas))
# # property_feed = mixture_property(T0, xi_gas=pd.Series(F0 / np.sum(F0), index=subs), Pt=P0)
# insula_sim = Insulation(Do, Din, 1, location=0)
# insula_sim.flux(Th=503, Pp=70, F_dict=gas, Tc=Tc)
# b = np.array([1,3,-1,-1,0])
# c = np.where(b>=0, b, 0)
# print(c)
comp = ["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]
const, cor = ChemicalConstantsPackage.from_IDs(comp)
print(const)


class ConstDB:
    def __init__(self, CASs=None, names=None, MWs=None, Tms=None, Tbs=None,
                 # Critical state points
                 Tcs=None, Pcs=None, Vcs=None, omegas=None,
                 Zcs=None, rhocs=None, rhocs_mass=None,
                 # Phase change enthalpy
                 Hfus_Tms=None, Hfus_Tms_mass=None, Hvap_Tbs=None,
                 Hvap_Tbs_mass=None,
                 # Standard values
                 Vml_STPs=None, rhol_STPs=None, rhol_STPs_mass=None,
                 Vml_60Fs=None, rhol_60Fs=None, rhol_60Fs_mass=None,
                 Vmg_STPs=None, rhog_STPs=None, rhog_STPs_mass=None,
                 # Reaction (ideal gas)
                 Hfgs=None, Hfgs_mass=None, Gfgs=None, Gfgs_mass=None,
                 Sfgs=None, Sfgs_mass=None, S0gs=None, S0gs_mass=None,
                 Hf_STPs=None, Hf_STPs_mass=None,

                 # Triple point
                 Tts=None, Pts=None, Hsub_Tts=None, Hsub_Tts_mass=None,
                 # Combustion
                 Hcs=None, Hcs_mass=None, Hcs_lower=None, Hcs_lower_mass=None,
                 # Fire safety
                 Tflashs=None, Tautoignitions=None, LFLs=None, UFLs=None,
                 # Other safety
                 TWAs=None, STELs=None, Ceilings=None, Skins=None,
                 Carcinogens=None, legal_statuses=None, economic_statuses=None,
                 # Environmental
                 GWPs=None, ODPs=None, logPs=None,
                 Psat_298s=None, Hvap_298s=None, Hvap_298s_mass=None,
                 Vml_Tms=None, rhos_Tms=None, Vms_Tms=None, rhos_Tms_mass=None,

                 # Analytical
                 sigma_STPs=None, sigma_Tbs=None, sigma_Tms=None,
                 RIs=None, RI_Ts=None, conductivities=None,
                 conductivity_Ts=None,
                 # Odd constants
                 charges=None, dipoles=None, Stockmayers=None,
                 molecular_diameters=None, Van_der_Waals_volumes=None,
                 Van_der_Waals_areas=None, Parachors=None, StielPolars=None,
                 atomss=None, atom_fractions=None,
                 similarity_variables=None, phase_STPs=None,
                 solubility_parameters=None,
                 # Other identifiers
                 PubChems=None, formulas=None, smiless=None, InChIs=None,
                 InChI_Keys=None, aliases=None,
                 # Groups
                 UNIFAC_groups=None, UNIFAC_Dortmund_groups=None,
                 PSRK_groups=None, UNIFAC_Rs=None, UNIFAC_Qs=None,
                 ):
        self.N = N = len(MWs)
        self.cmps = range(N)

        empty_list = [None] * N

        if atom_fractions is None: atom_fractions = empty_list
        if atomss is None: atomss = empty_list
        if Carcinogens is None: Carcinogens = empty_list
        if CASs is None: CASs = empty_list
        if aliases is None: aliases = empty_list
        if Ceilings is None: Ceilings = empty_list
        if charges is None: charges = empty_list
        if conductivities is None: conductivities = empty_list
        if dipoles is None: dipoles = empty_list
        if economic_statuses is None: economic_statuses = empty_list
        if formulas is None: formulas = empty_list
        if Gfgs is None: Gfgs = empty_list
        if Gfgs_mass is None: Gfgs_mass = empty_list
        if GWPs is None: GWPs = empty_list
        if Hcs is None: Hcs = empty_list
        if Hcs_lower is None: Hcs_lower = empty_list
        if Hcs_lower_mass is None: Hcs_lower_mass = empty_list
        if Hcs_mass is None: Hcs_mass = empty_list
        if Hfgs is None: Hfgs = empty_list
        if Hfgs_mass is None: Hfgs_mass = empty_list
        if Hfus_Tms is None: Hfus_Tms = empty_list
        if Hfus_Tms_mass is None: Hfus_Tms_mass = empty_list
        if Hsub_Tts is None: Hsub_Tts = empty_list
        if Hsub_Tts_mass is None: Hsub_Tts_mass = empty_list
        if Hvap_298s is None: Hvap_298s = empty_list
        if Hvap_298s_mass is None: Hvap_298s_mass = empty_list
        if Hvap_Tbs is None: Hvap_Tbs = empty_list
        if Hvap_Tbs_mass is None: Hvap_Tbs_mass = empty_list
        if InChI_Keys is None: InChI_Keys = empty_list
        if InChIs is None: InChIs = empty_list
        if legal_statuses is None: legal_statuses = empty_list
        if LFLs is None: LFLs = empty_list
        if logPs is None: logPs = empty_list
        if molecular_diameters is None: molecular_diameters = empty_list
        if names is None: names = empty_list
        if ODPs is None: ODPs = empty_list
        if omegas is None: omegas = empty_list
        if Parachors is None: Parachors = empty_list
        if Pcs is None: Pcs = empty_list
        if phase_STPs is None: phase_STPs = empty_list
        if Psat_298s is None: Psat_298s = empty_list
        if PSRK_groups is None: PSRK_groups = empty_list
        if Pts is None: Pts = empty_list
        if PubChems is None: PubChems = empty_list
        if rhocs is None: rhocs = empty_list
        if rhocs_mass is None: rhocs_mass = empty_list
        if rhol_STPs is None: rhol_STPs = empty_list
        if rhol_STPs_mass is None: rhol_STPs_mass = empty_list
        if RIs is None: RIs = empty_list
        if S0gs is None: S0gs = empty_list
        if S0gs_mass is None: S0gs_mass = empty_list
        if Sfgs is None: Sfgs = empty_list
        if Sfgs_mass is None: Sfgs_mass = empty_list
        if similarity_variables is None: similarity_variables = empty_list
        if Skins is None: Skins = empty_list
        if smiless is None: smiless = empty_list
        if STELs is None: STELs = empty_list
        if StielPolars is None: StielPolars = empty_list
        if Stockmayers is None: Stockmayers = empty_list
        if solubility_parameters is None: solubility_parameters = empty_list
        if Tautoignitions is None: Tautoignitions = empty_list
        if Tbs is None: Tbs = empty_list
        if Tcs is None: Tcs = empty_list
        if Tflashs is None: Tflashs = empty_list
        if Tms is None: Tms = empty_list
        if Tts is None: Tts = empty_list
        if TWAs is None: TWAs = empty_list
        if UFLs is None: UFLs = empty_list
        if UNIFAC_Dortmund_groups is None: UNIFAC_Dortmund_groups = empty_list
        if UNIFAC_groups is None: UNIFAC_groups = empty_list
        if UNIFAC_Rs is None: UNIFAC_Rs = empty_list
        if UNIFAC_Qs is None: UNIFAC_Qs = empty_list
        if Van_der_Waals_areas is None: Van_der_Waals_areas = empty_list
        if Van_der_Waals_volumes is None: Van_der_Waals_volumes = empty_list
        if Vcs is None: Vcs = empty_list
        if Vml_STPs is None: Vml_STPs = empty_list
        if Vml_Tms is None: Vml_Tms = empty_list
        if rhos_Tms is None: rhos_Tms = empty_list
        if Vms_Tms is None: Vms_Tms = empty_list
        if Zcs is None: Zcs = empty_list
        if Vml_60Fs is None: Vml_60Fs = empty_list
        if rhol_60Fs is None: rhol_60Fs = empty_list
        if rhol_60Fs_mass is None: rhol_60Fs_mass = empty_list
        if RI_Ts is None: RI_Ts = empty_list
        if conductivity_Ts is None: conductivity_Ts = empty_list
        if rhos_Tms_mass is None: rhos_Tms_mass = empty_list
        if Vmg_STPs is None: Vmg_STPs = empty_list
        if rhog_STPs is None: rhog_STPs = empty_list
        if rhog_STPs_mass is None: rhog_STPs_mass = empty_list
        if sigma_STPs is None: sigma_STPs = empty_list
        if sigma_Tbs is None: sigma_Tbs = empty_list
        if sigma_Tms is None: sigma_Tms = empty_list
        if Hf_STPs is None: Hf_STPs = empty_list
        if Hf_STPs_mass is None: Hf_STPs_mass = empty_list

        self.atom_fractions = atom_fractions
        self.atomss = atomss
        self.Carcinogens = Carcinogens
        self.CASs = CASs
        self.Ceilings = Ceilings
        self.charges = charges
        self.conductivities = conductivities
        self.dipoles = dipoles
        self.economic_statuses = economic_statuses
        self.formulas = formulas
        self.Gfgs = Gfgs
        self.Gfgs_mass = Gfgs_mass
        self.GWPs = GWPs
        self.Hcs = Hcs
        self.Hcs_lower = Hcs_lower
        self.Hcs_lower_mass = Hcs_lower_mass
        self.Hcs_mass = Hcs_mass
        self.Hfgs = Hfgs
        self.Hfgs_mass = Hfgs_mass
        self.Hfus_Tms = Hfus_Tms
        self.Hfus_Tms_mass = Hfus_Tms_mass
        self.Hsub_Tts = Hsub_Tts
        self.Hsub_Tts_mass = Hsub_Tts_mass
        self.Hvap_298s = Hvap_298s
        self.Hvap_298s_mass = Hvap_298s_mass
        self.Hvap_Tbs = Hvap_Tbs
        self.Hvap_Tbs_mass = Hvap_Tbs_mass
        self.InChI_Keys = InChI_Keys
        self.InChIs = InChIs
        self.legal_statuses = legal_statuses
        self.aliases = aliases
        self.LFLs = LFLs
        self.logPs = logPs
        self.molecular_diameters = molecular_diameters
        self.MWs = MWs
        self.names = names
        self.ODPs = ODPs
        self.omegas = omegas
        self.Parachors = Parachors
        self.Pcs = Pcs
        self.phase_STPs = phase_STPs
        self.Psat_298s = Psat_298s
        self.PSRK_groups = PSRK_groups
        self.Pts = Pts
        self.PubChems = PubChems
        self.rhocs = rhocs
        self.rhocs_mass = rhocs_mass
        self.rhol_STPs = rhol_STPs
        self.rhol_STPs_mass = rhol_STPs_mass
        self.RIs = RIs
        self.S0gs = S0gs
        self.S0gs_mass = S0gs_mass
        self.Sfgs = Sfgs
        self.Sfgs_mass = Sfgs_mass
        self.similarity_variables = similarity_variables
        self.solubility_parameters = solubility_parameters
        self.Skins = Skins
        self.smiless = smiless
        self.STELs = STELs
        self.StielPolars = StielPolars
        self.Stockmayers = Stockmayers
        self.Tautoignitions = Tautoignitions
        self.Tbs = Tbs
        self.Tcs = Tcs
        self.Tflashs = Tflashs
        self.Tms = Tms
        self.Tts = Tts
        self.TWAs = TWAs
        self.UFLs = UFLs
        self.UNIFAC_Dortmund_groups = UNIFAC_Dortmund_groups
        self.UNIFAC_groups = UNIFAC_groups
        self.UNIFAC_Rs = UNIFAC_Rs
        self.UNIFAC_Qs = UNIFAC_Qs
        self.Van_der_Waals_areas = Van_der_Waals_areas
        self.Van_der_Waals_volumes = Van_der_Waals_volumes
        self.Vcs = Vcs
        self.Vml_STPs = Vml_STPs
        self.Vml_Tms = Vml_Tms
        self.rhos_Tms = rhos_Tms
        self.Vms_Tms = Vms_Tms
        self.Zcs = Zcs
        self.Vml_60Fs = Vml_60Fs
        self.rhol_60Fs = rhol_60Fs
        self.rhol_60Fs_mass = rhol_60Fs_mass
        self.conductivity_Ts = conductivity_Ts
        self.RI_Ts = RI_Ts
        self.rhos_Tms_mass = rhos_Tms_mass
        self.Vmg_STPs = Vmg_STPs
        self.rhog_STPs = rhog_STPs
        self.rhog_STPs_mass = rhog_STPs_mass
        self.sigma_STPs = sigma_STPs
        self.sigma_Tbs = sigma_Tbs
        self.sigma_Tms = sigma_Tms
        self.Hf_STPs = Hf_STPs
        self.Hf_STPs_mass = Hf_STPs_mass

        try:
            self.n_atoms = [sum(i.values()) for i in atomss]
        except:
            self.n_atoms = None


const_data = dict(atom_fractions=[{'C': 0.3333333333333333, 'O': 0.6666666666666666},
                                  {'H': 1.0},
                                  {'C': 0.16666666666666666, 'H': 0.6666666666666666, 'O': 0.16666666666666666},
                                  {'H': 0.6666666666666666, 'O': 0.3333333333333333},
                                  {'C': 0.5, 'O': 0.5}],
                  atomss=[{'C': 1, 'O': 2}, {'H': 2}, {'C': 1, 'H': 4, 'O': 1}, {'H': 2, 'O': 1}, {'C': 1, 'O': 1}],
                  Carcinogens=[{'International Agency for Research on Cancer': 'Unlisted',
                          'National Toxicology Program 13th Report on Carcinogens': 'Unlisted'},
                         {'International Agency for Research on Cancer': 'Unlisted',
                          'National Toxicology Program 13th Report on Carcinogens': 'Unlisted'},
                         {'International Agency for Research on Cancer': 'Unlisted',
                          'National Toxicology Program 13th Report on Carcinogens': 'Unlisted'},
                         {'International Agency for Research on Cancer': 'Unlisted',
                          'National Toxicology Program 13th Report on Carcinogens': 'Unlisted'},
                         {'International Agency for Research on Cancer': 'Unlisted',
                          'National Toxicology Program 13th Report on Carcinogens': 'Unlisted'}],
                  CASs=['124-38-9', '1333-74-0', '67-56-1', '7732-18-5', '630-08-0'],
                  charges=[0, 0, 0, 0, 0], conductivities=[None, None, 4.4e-05, 4e-06, None],
                  dipoles=[0.0, 0.0, 1.7, 1.85, 0.11], formulas=['CO2', 'H2', 'CH4O', 'H2O', 'CO'],
                  Gfgs=[-394338.635, 0.0, -162000.13, -228554.325, -137179.61],
                  Gfgs_mass=[-8960307.092786783, 0.0, -5055890.325967344, -12686692.9073542, -4897505.185629469],
                  GWPs=[1.0, None, None, None, None], Hcs=[0.0, -285825.0, -726724.0, 0.0, -282949.0],
                  Hcs_lower=[0.0, -241813.50400000002, -638701.008, 44011.496, -282949.0],
                  Hcs_lower_mass=[0.0, -119954314.74095681, -19933331.211109467, 2443009.2676883177, -10101677.609148128],
                  Hcs_mass=[0.0, -141786713.49485087, -22680456.12832713, 0.0, -10101677.609148128],
                  Hfgs=[-393474.0, 0.0, -200700.0, -241822.0, -110525.0],
                  Hfgs_mass=[-8940660.539201763, 0.0, -6263681.321870828, -13423160.783512661, -3945898.0867615608],
                  Hfus_Tms=[9020.0, 120.0, 3215.0, 6010.0, 833.0],
                  Hfus_Tms_mass=[204955.74819073154, 59527.35281862015, 100337.49601302797, 333605.69472136983,
                           29739.272619519386],
                  Hsub_Tts=[24435.40534684326, 1033.3903533945663, 45351.41225231173, 51065.16013858828, 7330.542813643393],
                  Hsub_Tts_mass=[555230.2422623129, 512624.93471564096, 1415380.1387407514, 2834547.125472836,
                           261710.69769987944],
                  Hvap_298s=[5265.543624363681, 0.0, 37457.19357385385, 43987.4460555689, 0.0],
                  Hvap_298s_mass=[119645.61343263797, 0.0, 1169008.0904745809, 2441674.2929096245, 0.0],
                  Hvap_Tbs=[16848.54074668973, 904.5470282563128, 35280.6908513165, 40650.94453927768, 6013.289911773105],
                  Hvap_Tbs_mass=[382838.7222460998, 448710.7507670659, 1101081.237210215, 2256470.315159003,
                           214682.9147976303],
                  InChI_Keys=['CURLTUGMZLYLDI-UHFFFAOYSA-N', 'UFHFLCQGNIYNRP-UHFFFAOYSA-N', 'OKKJLVBELUTLKV-UHFFFAOYSA-N',
                        'XLYOFNOQVPJJNP-UHFFFAOYSA-N', 'UGFAIRIUMAVXCW-UHFFFAOYSA-N'],
                  InChIs=['CO2/c2-1-3', 'H2/h1H', 'CH4O/c1-2/h2H,1H3', 'H2O/h1H2', 'CO/c1-2'],
                  LFLs=[None, 0.04, 0.06, None, 0.109],
                  logPs=[0.83, None, -0.74, -1.38, None],
                  molecular_diameters=[3.26192, 5.94111, 3.79957, 3.24681, 3.53562],
                  MWs=[44.0095, 2.01588, 32.04186, 18.01528, 28.0101],
                  names=['carbon dioxide', 'hydrogen', 'methanol', 'water', 'carbon monoxide'],
                  aliases=['CO2', 'H2', 'Methanol', 'H2O', 'carbon monoxide'],
                  omegas=[0.22394, -0.219, 0.5625, 0.3443, 0.0497],
                  Parachors=[7.372678555293655e-06, 0.0, 1.5746204696062765e-05, 9.365630122941035e-06, 0.0],
                  Pcs=[7377300.0, 1296400.0, 8215850.0, 22064000.0, 3494000.0], phase_STPs=['g', 'g', 'l', 'l', 'g'],
                  Psat_298s=[6434243.259502665, 6.684819223510374e+16, 16981.49712384057, 3169.9293388738784,
                       123835973922.04547],
                  PSRK_groups=[{}, {}, {15: 1}, {16: 1}, {}],
                  Pts=[517964.343354, 7357.81734541, 0.186349762324, 611.654771008, 15536.8538588],
                  PubChems=[280, 783, 887, 962, 281],
                  rhocs=[10624.9062999959, 15508.00000001163, 8785.16999996665, 17873.727995602905, 10850.000000001357],
                  rhocs_mass=[467.5968138096697, 31.262267040023453, 281.4931872151314, 322.00021448462513, 303.909585000038],
                  rhol_STPs=[20995.90611365801, 27287.30743919876, 24540.340077095232, 55344.57560175075, 20641.754273261144],
                  rhol_STPs_mass=[924.0193301090322, 55.007937320532, 786.3181411026746, 997.0480259467083, 578.177601369472],
                  RIs=[None, None, 1.3288, 1.3945, None], S0gs=[213.8, 130.7, 239.9, 188.8, 197.7],
                  S0gs_mass=[4858.042013656142, 64835.20844494711, 7487.080962216301, 10479.992539666328, 7058.16830357621],
                  Sfgs=[2.9000000000000314, 0.0, -129.8, -44.499999999999964, 89.39999999999996],
                  Sfgs_mass=[65.89486360899423, 0.0, -4050.9508499194494, -2470.1253602497413, 3191.705848961623],
                  similarity_variables=[0.0681671002851657, 0.9921225469770025, 0.18725504699165404, 0.16652530518537598,
                                  0.07140281541301173],
                  Skins=[False, False, True, None, False],
                  smiless=['C(=O)=O', '[HH]', 'CO', 'O', '[C-]#[O+]'],
                  STELs=[(30000.0, 'ppm'), None, (250.0, 'ppm'), None, None],
                  StielPolars=[0.0143868058231571, 0.009638368562742805, 0.03852800834983805, 0.02350430402645798,
                         0.001354442464295058],
                  Stockmayers=[500.71, 3.4503, 685.96, 501.01, 102.86],
                  Tautoignitions=[None, 833.15, 713.15, None, 880.15],
                  Tbs=[194.67, 20.3689085101, 337.632383296, 373.124295848, 81.6381829183],
                  Tcs=[304.1282, 33.145, 513.38, 647.096, 132.86], Tflashs=[None, None, 282.15, None, None],
                  Tms=[216.65, 13.99, 175.15, 273.15, 68.15], Tts=[216.592, 13.957, 175.61, 273.16, 68.16],
                  TWAs=[(5000.0, 'ppm'), None, (200.0, 'ppm'), None, (25.0, 'ppm')], UFLs=[None, 0.77, 0.36, None, 0.74],
                  UNIFAC_Dortmund_groups=[{}, {}, {15: 1}, {16: 1}, {}],
                  UNIFAC_groups=[{}, {}, {15: 1}, {16: 1}, {}],
                  Van_der_Waals_areas=[0.0, 0.0, 358000.0, 350000.0, 0.0],
                  Van_der_Waals_volumes=[0.0, 0.0, 2.1709787e-05, 1.39564e-05, 0.0],
                  Vcs=[9.41184770731e-05, 6.44828475625e-05, 0.000113828190007, 5.59480372671e-05, 9.21658986175e-05],
                  Vml_STPs=[4.762833261811415e-05, 3.664707491672411e-05, 4.074923154522018e-05, 1.806861809178579e-05,
                      4.844549483351699e-05],
                  Vml_Tms=[3.7351689383750086e-05, 2.616094079172244e-05, 3.542096205028121e-05, 1.8018099870410792e-05,
                     3.2965274497311146e-05],
                  Zcs=[0.27458794013625015, 0.3033409353716562, 0.21909335258418086, 0.22943845208106295, 0.2915175660594411],
                  UNIFAC_Rs=[0.0, 0.0, 1.4311, 0.92, 0.0], UNIFAC_Qs=[0.0, 0.0, 1.432, 1.4, 0.0],
                  rhos_Tms=[29984.109811216116, 42797.46022768705, 31631.526804340934, 62160.04270307297, 33975.857664093164],
                  Vms_Tms=[3.335099845538624e-05, 2.3365872523273424e-05, 3.1614028819587856e-05, 1.6087505035619347e-05,
                     2.9432663919381608e-05],
                  rhos_Tms_mass=[1319.5856807367159, 86.27454412378978, 1013.5329534509395, 1119.8305741078166,
                           951.667170757016],
                  solubility_parameters=[7648.981011950944, None, 29298.08560457511, 47929.841529177254, None],
                  Vml_60Fs=[4.762833261811415e-05, 3.664707491672411e-05, 4.029593850336912e-05, 1.8032991611014203e-05,
                      4.844549483351699e-05],
                  rhol_60Fs=[20995.90611365801, 27287.30743919876, 24816.396816676464, 55453.91588765667, 20641.754273261144],
                  rhol_60Fs_mass=[924.0193301090322, 55.007937320532, 795.1635125043929, 999.0178218125835, 578.177601369472],
                  conductivity_Ts=[None, None, 291.15, 291.15, None], RI_Ts=[None, None, 293.15, None, None],
                  Vmg_STPs=[0.024465403697038125, 0.024465403697038125, 0.024465403697038125, 0.024465403697038125,
                      0.024465403697038125],
                  rhog_STPs=[40.87404452357612, 40.87404452357612, 40.87404452357612, 40.87404452357612, 40.87404452357612],
                  rhog_STPs_mass=[1.7988462624603232, 0.08239716887418663, 1.3096804122581924, 0.7363573568246904,
                            1.1448860745098195],
                  sigma_STPs=[0.0005697108840068718, 0.0, 0.02214777216932284, 0.07197220523022962, 0.0],
                  sigma_Tms=[0.016480773148367444, 0.003026548437349704, 0.035137430555304106, 0.07564766822989494,
                       0.012448454316005263],
                  sigma_Tbs=[0.021829936817407598, 0.0019116526036083092, 0.0188130778659449, 0.058916822384573755,
                       0.009518632965656058],
                  Hf_STPs=[-393474.0, 0.0, -238400.0, -285825.0, -110525.0],
                  Hf_STPs_mass=[-8940660.539201763, 0.0, -7440267.200468387, -15865698.451536695, -3945898.0867615608])

a= ConstDB(**const_data)
print(a.Sfgs)